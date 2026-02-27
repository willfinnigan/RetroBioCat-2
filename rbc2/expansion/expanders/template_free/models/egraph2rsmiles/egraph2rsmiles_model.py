"""EGraph2RSMILES model and inference wrapper.

Ported from https://gitlab.com/mlpds_mit/askcosv2/retro/egraph2rsmiles
"""

import logging
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
from rbc2.expansion.expanders.template_free.utils.chem_utils import ATOM_FDIMS, BOND_FDIMS
from rbc2.expansion.expanders.template_free.utils.data_utils import (
    EG2RSBatch, create_batch_from_smi
)
from rbc2.configs.data_path import path_to_data_folder

logger = logging.getLogger(__name__)


class EGraph2RSMILES(nn.Module):
    """EGraph2RSMILES neural network model."""

    def __init__(self, args, vocab: dict):
        super().__init__()
        self.args = args
        self.vocab = vocab
        self.vocab_size = len(self.vocab)
        self.n_latent = args.n_latent

        while args.enable_amp and not self.vocab_size % 8 == 0:
            self.vocab_size += 1

        # pad + node_type + ATOM_FDIM + BOND_FDIM
        self.feature_lookup = nn.Embedding(
            1 + 2 + sum(ATOM_FDIMS) + sum(BOND_FDIMS),
            args.encoder_hidden_size,
            padding_idx=0)
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
            attention_dropout=args.dropout,
            embeddings=self.decoder_embeddings,
            max_relative_positions=args.max_relative_positions,
            full_context_alignment=False,
            alignment_layer=-3,
            alignment_heads=0,
            latent_embedding=self.latent_embedding,
            shared_decoder_layer=args.shared_decoder_layer)

        if not args.encoder_hidden_size == args.decoder_hidden_size:
            self.bridge_layer = nn.Linear(args.encoder_hidden_size, args.decoder_hidden_size, bias=True)

        self.output_layer = nn.Linear(args.decoder_hidden_size, self.vocab_size, bias=True)

        self.criterion = nn.CrossEntropyLoss(
            ignore_index=self.vocab["_PAD"],
            reduction="none")

    def encode_and_reshape(self, reaction_batch: EG2RSBatch):
        memory_lengths = reaction_batch.node_lengths
        max_length = max(memory_lengths)

        padded_memory_bank = self.feature_lookup(reaction_batch.fnode)
        padded_memory_bank = padded_memory_bank.sum(dim=2)
        padded_memory_bank = padded_memory_bank.transpose(0, 1).contiguous()
        padded_memory_bank = self.attention_encoder(
            padded_memory_bank, memory_lengths, reaction_batch.distances)

        if not self.args.encoder_hidden_size == self.args.decoder_hidden_size:
            padded_memory_bank = self.bridge_layer(padded_memory_bank)

        self.decoder.state["src"] = np.zeros(max_length)

        return padded_memory_bank, memory_lengths


@dataclass
class EGraph2RSmilesConfig:
    """Configuration for EGraph2RSMILES model architecture."""
    embed_size: int = 256
    shared_encoder_layer: int = 0
    shared_decoder_layer: int = 0
    encoder_num_layers: int = 4
    encoder_hidden_size: int = 256
    encoder_attn_heads: int = 8
    encoder_filter_size: int = 2048
    decoder_num_layers: int = 4
    decoder_hidden_size: int = 256
    decoder_attn_heads: int = 8
    decoder_filter_size: int = 2048
    rel_pos_buckets: int = 31
    max_relative_positions: int = 0
    n_latent: int = 0
    dropout: float = 0.0
    enable_amp: bool = False


class EGraph2RSmilesModel(GraphInferenceModel):
    """Inference wrapper for EGraph2RSMILES."""

    def __init__(self,
                 checkpoint_path: Optional[str] = None,
                 vocab_path: Optional[str] = None,
                 beam_size: int = 10,
                 max_length: int = 512,
                 temperature: float = 1.0):

        data_dir = os.path.join(path_to_data_folder, 'template_free', 'egraph2rsmiles')
        super().__init__(
            checkpoint_path=checkpoint_path or os.path.join(data_dir, 'model.pt'),
            vocab_path=vocab_path or os.path.join(data_dir, 'vocab.txt'),
            beam_size=beam_size,
            max_length=max_length,
            temperature=temperature)
        self.rel_pos_buckets = EGraph2RSmilesConfig.rel_pos_buckets

    def _ensure_weights_exist(self):
        from rbc2.configs.download_data_files.download_template_free import ensure_egraph2rsmiles_exists
        ensure_egraph2rsmiles_exists()

    def _create_config(self):
        return EGraph2RSmilesConfig()

    def _create_model(self, config, vocab):
        return EGraph2RSMILES(config, vocab)

    def _create_batch(self, smi):
        return create_batch_from_smi(
            smi, rel_pos_buckets=self.rel_pos_buckets, device=self._device)

    def _normalize_score(self, beam_score):
        # Linear mapping from log-prob: score of 0.0 -> 1.0, score of -5.0 -> 0.0
        return max(0.0, min(1.0, (beam_score + 5.0) / 5.0))

    def _on_config_loaded(self, config):
        self.rel_pos_buckets = config.rel_pos_buckets
