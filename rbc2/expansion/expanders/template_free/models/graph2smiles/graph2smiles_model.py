"""Graph2SMILES model and inference wrapper.

Ported from https://gitlab.com/mlpds_mit/askcosv2/retro/graph2smiles
"""

import logging
import math
import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import torch
import torch.nn as nn

from rbc2.expansion.expanders.template_free.models.graph_inference_base import GraphInferenceModel
from rbc2.expansion.expanders.template_free.models.model_interface import TemplateFreeResult
from rbc2.expansion.expanders.template_free.models.shared.attention_xl import AttnEncoderXL
from rbc2.expansion.expanders.template_free.models.shared.transformer_decoder import TransformerDecoder
from rbc2.expansion.expanders.template_free.models.shared.embeddings import Embeddings
from rbc2.expansion.expanders.template_free.models.graph2smiles.graphfeat import GraphFeatEncoder
from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
    ATOM_FDIM, BOND_FDIM, create_g2s_batch_from_smi
)
from rbc2.configs.data_path import path_to_data_folder

logger = logging.getLogger(__name__)


# --- Model config ---

@dataclass
class Graph2SmilesConfig:
    """Configuration for Graph2SMILES model architecture."""
    embed_size: int = 256
    encoder_num_layers: int = 4
    encoder_hidden_size: int = 256
    encoder_attn_heads: int = 8
    encoder_filter_size: int = 2048
    attn_enc_num_layers: int = 4
    attn_enc_hidden_size: int = 256
    attn_enc_heads: int = 8
    attn_enc_filter_size: int = 2048
    decoder_num_layers: int = 4
    decoder_hidden_size: int = 256
    decoder_attn_heads: int = 8
    decoder_filter_size: int = 2048
    rel_pos_buckets: int = 31
    rel_pos: str = 'emb_only'
    max_relative_positions: int = 0
    n_latent: int = 0
    dropout: float = 0.0
    attn_dropout: float = 0.0
    enable_amp: bool = False
    mpn_type: str = 'dgcn'
    dgat_attn_heads: int = 4
    shared_encoder_layer: int = 0
    shared_decoder_layer: int = 0
    encoder_positional_encoding: str = 'none'
    encoder_emb_scale: str = 'sqrt'
    compute_graph_distance: bool = True
    mask_rel_chirality: int = 0


# --- Graph2SMILES neural network ---

class Graph2SMILES(nn.Module):
    """Graph2SMILES neural network model."""

    def __init__(self, args, vocab: dict):
        super().__init__()
        self.args = args
        self.vocab = vocab
        self.vocab_size = len(self.vocab)
        self.n_latent = int(args.n_latent)

        while args.enable_amp and not self.vocab_size % 8 == 0:
            self.vocab_size += 1

        self.encoder = GraphFeatEncoder(
            args,
            n_atom_feat=sum(ATOM_FDIM),
            n_bond_feat=BOND_FDIM)

        self.attention_encoder = None
        if args.attn_enc_num_layers > 0:
            self.attention_encoder = AttnEncoderXL(args)

        self.latent_embedding = None
        if self.n_latent > 1:
            self.latent_embedding = nn.Embedding(self.n_latent, args.decoder_hidden_size)

        self.decoder_embeddings = Embeddings(
            word_vec_size=args.embed_size,
            word_vocab_size=self.vocab_size,
            word_padding_idx=self.vocab["_PAD"],
            position_encoding=True,
            dropout=args.dropout)

        self.decoder = TransformerDecoder(
            num_layers=args.decoder_num_layers,
            d_model=args.decoder_hidden_size,
            heads=args.decoder_attn_heads,
            d_ff=args.decoder_filter_size,
            copy_attn=False,
            self_attn_type="scaled-dot",
            dropout=args.dropout,
            attention_dropout=getattr(args, 'attn_dropout', args.dropout),
            embeddings=self.decoder_embeddings,
            max_relative_positions=args.max_relative_positions,
            full_context_alignment=False,
            alignment_layer=-3,
            alignment_heads=0,
            latent_embedding=self.latent_embedding,
            shared_decoder_layer=getattr(args, 'shared_decoder_layer', 0))

        if not args.attn_enc_hidden_size == args.decoder_hidden_size:
            self.bridge_layer = nn.Linear(args.attn_enc_hidden_size, args.decoder_hidden_size, bias=True)

        self.output_layer = nn.Linear(args.decoder_hidden_size, self.vocab_size, bias=True)

        self.criterion = nn.CrossEntropyLoss(
            ignore_index=self.vocab["_PAD"],
            reduction="none")

    def encode_and_reshape(self, reaction_batch):
        hatom, _ = self.encoder(reaction_batch)

        if not self.args.attn_enc_hidden_size == self.args.decoder_hidden_size:
            hatom = self.bridge_layer(hatom)

        atom_scope = reaction_batch.atom_scope
        memory_lengths = [scope[-1][0] + scope[-1][1] - scope[0][0]
                          for scope in atom_scope]

        memory_bank = torch.split(hatom, [1] + memory_lengths, dim=0)
        padded_memory_bank = []
        max_length = max(memory_lengths)

        for length, h in zip(memory_lengths, memory_bank[1:]):
            m = nn.ZeroPad2d((0, 0, 0, max_length - length))
            padded_memory_bank.append(m(h))

        padded_memory_bank = torch.stack(padded_memory_bank, dim=1)

        memory_lengths = torch.tensor(memory_lengths, dtype=torch.long,
                                       device=padded_memory_bank.device)

        if self.attention_encoder is not None:
            padded_memory_bank = self.attention_encoder(
                padded_memory_bank, memory_lengths, reaction_batch.distances)

        self.decoder.state["src"] = np.zeros(max_length)

        return padded_memory_bank, memory_lengths


# --- Inference wrapper ---

class Graph2SmilesModel(GraphInferenceModel):
    """Inference wrapper for Graph2SMILES."""

    def __init__(self,
                 checkpoint_path: Optional[str] = None,
                 vocab_path: Optional[str] = None,
                 beam_size: int = 10,
                 max_length: int = 512,
                 temperature: float = 1.0):

        data_dir = os.path.join(path_to_data_folder, 'template_free', 'graph2smiles')
        super().__init__(
            checkpoint_path=checkpoint_path or os.path.join(data_dir, 'model.pt'),
            vocab_path=vocab_path or os.path.join(data_dir, 'vocab.txt'),
            beam_size=beam_size,
            max_length=max_length,
            temperature=temperature)
        self.rel_pos_buckets = Graph2SmilesConfig.rel_pos_buckets

    def _ensure_weights_exist(self):
        from rbc2.configs.download_data_files.download_template_free import ensure_graph2smiles_exists
        ensure_graph2smiles_exists()

    def _create_config(self):
        return Graph2SmilesConfig()

    def _create_model(self, config, vocab):
        return Graph2SMILES(config, vocab)

    def _create_batch(self, smi):
        return create_g2s_batch_from_smi(
            smi, rel_pos_buckets=self.rel_pos_buckets,
            compute_graph_distance=True, device=self._device)

    def _on_config_loaded(self, config):
        self.rel_pos_buckets = config.rel_pos_buckets
