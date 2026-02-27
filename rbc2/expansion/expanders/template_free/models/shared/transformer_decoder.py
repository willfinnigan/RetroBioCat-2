"""Adapted from onmt.decoder.transformer"""

import numpy as np
import torch
import torch.nn as nn
from rbc2.expansion.expanders.template_free.models.shared.attention_onmt import MultiHeadedAttention
from rbc2.expansion.expanders.template_free.models.shared.position_ffn import PositionwiseFeedForward
from rbc2.expansion.expanders.template_free.models.shared.misc import sequence_mask


class TransformerDecoderLayer(nn.Module):
    """Transformer Decoder layer block in Pre-Norm style."""

    def __init__(self, d_model, heads, d_ff, dropout, attention_dropout,
                 self_attn_type="scaled-dot", max_relative_positions=0,
                 full_context_alignment=False, alignment_heads=0):
        super().__init__()

        assert self_attn_type == "scaled-dot"
        self.self_attn = MultiHeadedAttention(
            heads, d_model, dropout=attention_dropout,
            max_relative_positions=max_relative_positions)
        self.context_attn = MultiHeadedAttention(
            heads, d_model, dropout=attention_dropout)
        self.feed_forward = PositionwiseFeedForward(d_model, d_ff, dropout)
        self.layer_norm_1 = nn.LayerNorm(d_model, eps=1e-6)
        self.layer_norm_2 = nn.LayerNorm(d_model, eps=1e-6)
        self.drop = nn.Dropout(dropout)
        self.full_context_alignment = full_context_alignment
        self.alignment_heads = alignment_heads

    def forward(self, *args, **kwargs):
        with_align = kwargs.pop('with_align', False)
        output, attns = self._forward(*args, **kwargs)
        top_attn = attns[:, 0, :, :].contiguous()
        attn_align = None
        if with_align:
            if self.full_context_alignment:
                _, attns = self._forward(*args, **kwargs, future=True)
            if self.alignment_heads > 0:
                attns = attns[:, :self.alignment_heads, :, :].contiguous()
            attn_align = attns.mean(dim=1)
        return output, top_attn, attn_align

    def _forward(self, inputs, memory_bank, src_pad_mask, tgt_pad_mask,
                 layer_cache=None, step=None, future=False):
        dec_mask = None

        if step is None:
            tgt_len = tgt_pad_mask.size(-1)
            if not future:
                future_mask = torch.ones(
                    [tgt_len, tgt_len],
                    device=tgt_pad_mask.device,
                    dtype=torch.uint8)
                future_mask = future_mask.triu_(1).view(1, tgt_len, tgt_len)
                try:
                    future_mask = future_mask.bool()
                except AttributeError:
                    pass
                dec_mask = torch.gt(tgt_pad_mask + future_mask, 0)
            else:
                dec_mask = tgt_pad_mask

        input_norm = self.layer_norm_1(inputs)
        query, _ = self.self_attn(input_norm, input_norm, input_norm,
                                  mask=dec_mask,
                                  layer_cache=layer_cache,
                                  attn_type="self")
        query = self.drop(query) + inputs

        query_norm = self.layer_norm_2(query)
        mid, attns = self.context_attn(memory_bank, memory_bank, query_norm,
                                       mask=src_pad_mask,
                                       layer_cache=layer_cache,
                                       attn_type="context")
        output = self.feed_forward(self.drop(mid) + query)
        return output, attns

    def update_dropout(self, dropout, attention_dropout):
        self.self_attn.update_dropout(attention_dropout)
        self.context_attn.update_dropout(attention_dropout)
        self.feed_forward.update_dropout(dropout)
        self.drop.p = dropout


class TransformerDecoder(nn.Module):
    """The Transformer decoder from "Attention is All You Need"."""

    def __init__(self, num_layers, d_model, heads, d_ff,
                 copy_attn, self_attn_type, dropout, attention_dropout,
                 embeddings, max_relative_positions,
                 full_context_alignment, alignment_layer,
                 alignment_heads, latent_embedding, shared_decoder_layer):
        super(TransformerDecoder, self).__init__()

        self.embeddings = embeddings
        self.latent_embedding = latent_embedding
        self.num_layers = num_layers
        self.shared_decoder_layer = shared_decoder_layer

        self.state = {}

        if shared_decoder_layer:
            self.transformer_layer = TransformerDecoderLayer(
                d_model, heads, d_ff, dropout,
                attention_dropout, self_attn_type=self_attn_type,
                max_relative_positions=max_relative_positions,
                full_context_alignment=full_context_alignment,
                alignment_heads=alignment_heads)
        else:
            self.transformer_layers = nn.ModuleList(
                [TransformerDecoderLayer(
                    d_model, heads, d_ff, dropout,
                    attention_dropout, self_attn_type=self_attn_type,
                    max_relative_positions=max_relative_positions,
                    full_context_alignment=full_context_alignment,
                    alignment_heads=alignment_heads)
                 for i in range(num_layers)])

        self._copy = copy_attn
        self.layer_norm = nn.LayerNorm(d_model, eps=1e-6)
        self.alignment_layer = alignment_layer

    def init_state(self, src, memory_bank, enc_hidden):
        self.state["src"] = src
        self.state["cache"] = None

    def map_state(self, fn):
        def _recursive_map(struct, batch_dim=0):
            for k, v in struct.items():
                if v is not None:
                    if isinstance(v, dict):
                        _recursive_map(v)
                    else:
                        struct[k] = fn(v, batch_dim)

        self.state["src"] = fn(self.state["src"], 1)
        if self.state["cache"] is not None:
            _recursive_map(self.state["cache"])

    def detach_state(self):
        self.state["src"] = self.state["src"].detach()

    def forward(self, tgt, memory_bank, step=None, **kwargs):
        """Decode, possibly stepwise."""
        if step == 0:
            self._init_cache(memory_bank)

        tgt_words = tgt[:, :, 0].transpose(0, 1)

        emb = self.embeddings(tgt, step=step)
        if self.latent_embedding is not None:
            latent_emb = self.latent_embedding(kwargs["latent_classes"])
            if kwargs["mode"] == "train":
                latent_emb = torch.unsqueeze(latent_emb, dim=0)
            elif kwargs["mode"] == "predict":
                latent_emb = torch.reshape(latent_emb, [1, 1, -1])
            else:
                raise NotImplementedError
            emb = emb + latent_emb
        assert emb.dim() == 3

        output = emb.transpose(0, 1).contiguous()
        src_memory_bank = memory_bank.transpose(0, 1).contiguous()

        pad_idx = self.embeddings.word_padding_idx
        src_lens = kwargs["memory_lengths"]
        src_max_len = self.state["src"].shape[0]
        src_pad_mask = ~sequence_mask(src_lens, src_max_len).unsqueeze(1)
        tgt_pad_mask = tgt_words.data.eq(pad_idx).unsqueeze(1)

        with_align = kwargs.pop('with_align', False)
        attn_aligns = []

        if self.shared_decoder_layer:
            for i in range(self.num_layers):
                layer_cache = self.state["cache"]["layer_{}".format(i)] \
                    if step is not None else None
                output, attn, attn_align = self.transformer_layer(
                    output, src_memory_bank, src_pad_mask, tgt_pad_mask,
                    layer_cache=layer_cache, step=step, with_align=with_align)
                if attn_align is not None:
                    attn_aligns.append(attn_align)
        else:
            for i, layer in enumerate(self.transformer_layers):
                layer_cache = self.state["cache"]["layer_{}".format(i)] \
                    if step is not None else None
                output, attn, attn_align = layer(
                    output, src_memory_bank, src_pad_mask, tgt_pad_mask,
                    layer_cache=layer_cache, step=step, with_align=with_align)
                if attn_align is not None:
                    attn_aligns.append(attn_align)

        output = self.layer_norm(output)
        dec_outs = output.transpose(0, 1).contiguous()
        attn = attn.transpose(0, 1).contiguous()

        attns = {"std": attn}
        if self._copy:
            attns["copy"] = attn
        if with_align:
            attns["align"] = attn_aligns[self.alignment_layer]

        return dec_outs, attns

    def _init_cache(self, memory_bank):
        self.state["cache"] = {}
        batch_size = memory_bank.size(1)
        depth = memory_bank.size(-1)

        num_layers = self.num_layers
        for i in range(num_layers):
            layer_cache = {"memory_keys": None, "memory_values": None}
            layer_cache["self_keys"] = None
            layer_cache["self_values"] = None
            self.state["cache"]["layer_{}".format(i)] = layer_cache

    def update_dropout(self, dropout, attention_dropout):
        self.embeddings.update_dropout(dropout)
        if self.shared_decoder_layer:
            self.transformer_layer.update_dropout(dropout, attention_dropout)
        else:
            for layer in self.transformer_layers:
                layer.update_dropout(dropout, attention_dropout)
