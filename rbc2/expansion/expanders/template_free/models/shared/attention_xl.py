import math
import torch
import torch.nn as nn
from rbc2.expansion.expanders.template_free.models.shared.positional_encoding import PositionalEncoding
from rbc2.expansion.expanders.template_free.models.shared.position_ffn import PositionwiseFeedForward
from rbc2.expansion.expanders.template_free.models.shared.misc import sequence_mask


class MultiHeadedRelAttention(nn.Module):
    def __init__(self, args, head_count, model_dim, dropout, rel_pos_buckets, u, v):
        super().__init__()
        self.args = args

        assert model_dim % head_count == 0
        self.dim_per_head = model_dim // head_count
        self.model_dim = model_dim
        self.head_count = head_count

        self.linear_keys = nn.Linear(model_dim, model_dim)
        self.linear_values = nn.Linear(model_dim, model_dim)
        self.linear_query = nn.Linear(model_dim, model_dim)

        self.softmax = nn.Softmax(dim=-1)
        self.dropout = nn.Dropout(dropout)
        self.final_linear = nn.Linear(model_dim, model_dim)

        self.rel_pos_buckets = rel_pos_buckets

        rel_pos = getattr(args, 'rel_pos', 'emb_only')
        if rel_pos == "emb_only":
            self.relative_pe = nn.Embedding(
                rel_pos_buckets + 1,
                model_dim,
                padding_idx=rel_pos_buckets
            )
        elif rel_pos == "enc_only":
            pe = _get_sin_encodings(rel_pos_buckets, model_dim)
            self.relative_pe = nn.Embedding.from_pretrained(
                embeddings=pe, freeze=True, padding_idx=rel_pos_buckets
            )
        else:
            self.relative_pe = None

        self.u = u
        self.v = v

    def forward(self, inputs, mask, distances):
        batch_size = inputs.size(0)
        dim_per_head = self.dim_per_head
        head_count = self.head_count

        def shape(x):
            return x.view(batch_size, -1, head_count, dim_per_head).transpose(1, 2)

        def unshape(x):
            return x.transpose(1, 2).contiguous().view(batch_size, -1, head_count * dim_per_head)

        query = self.linear_query(inputs)
        key = self.linear_keys(inputs)
        value = self.linear_values(inputs)

        key = shape(key)
        value = shape(value)
        query = shape(query)

        key_len = key.size(2)
        query_len = query.size(2)

        query = query / math.sqrt(dim_per_head)

        if self.relative_pe is None:
            scores = torch.matmul(query, key.transpose(2, 3))
        else:
            u = self.u.reshape(1, head_count, 1, dim_per_head)
            a_c = torch.matmul(query + u, key.transpose(2, 3))

            rel_emb = self.relative_pe(distances)
            rel_emb = rel_emb.reshape(
                batch_size, query_len, key_len, head_count, dim_per_head)

            query = query.unsqueeze(-2)
            rel_emb_t = rel_emb.permute(0, 3, 1, 4, 2)

            v = self.v.reshape(1, head_count, 1, 1, dim_per_head)
            b_d = torch.matmul(query + v, rel_emb_t).squeeze(-2)

            scores = a_c + b_d

        scores = scores.float()
        mask = mask.unsqueeze(1)
        scores = scores.masked_fill(mask, -1e18)

        attn = self.softmax(scores)
        drop_attn = self.dropout(attn)

        context_original = torch.matmul(drop_attn, value)
        context = unshape(context_original)

        output = self.final_linear(context)
        attns = attn.view(batch_size, head_count, query_len, key_len)

        return output, attns


def _get_sin_encodings(rel_pos_buckets, model_dim):
    pe = torch.zeros(rel_pos_buckets + 1, model_dim)
    position = torch.arange(0, rel_pos_buckets).unsqueeze(1)
    div_term = torch.exp((torch.arange(0, model_dim, 2, dtype=torch.float) *
                          -(math.log(10000.0) / model_dim)))
    pe[:-1, 0::2] = torch.sin(position.float() * div_term)
    pe[:-1, 1::2] = torch.cos(position.float() * div_term)
    return pe


class SALayerXL(nn.Module):
    """A single layer of the self-attention encoder."""

    def __init__(self, args, d_model, heads, d_ff, dropout, attention_dropout, rel_pos_buckets, u, v):
        super().__init__()
        self.self_attn = MultiHeadedRelAttention(
            args, heads, d_model, dropout=attention_dropout,
            rel_pos_buckets=rel_pos_buckets, u=u, v=v)
        self.feed_forward = PositionwiseFeedForward(d_model, d_ff, dropout)
        self.layer_norm = nn.LayerNorm(d_model, eps=1e-6)
        self.dropout = nn.Dropout(dropout)

    def forward(self, inputs, mask, distances):
        input_norm = self.layer_norm(inputs)
        context, _ = self.self_attn(input_norm, mask=mask, distances=distances)
        out = self.dropout(context) + inputs
        return self.feed_forward(out)


class AttnEncoderXL(nn.Module):
    def __init__(self, args):
        super().__init__()
        self.args = args

        self.num_layers = getattr(args, 'attn_enc_num_layers', args.encoder_num_layers)
        self.d_model = getattr(args, 'attn_enc_hidden_size', args.encoder_hidden_size)
        self.heads = getattr(args, 'attn_enc_heads', args.encoder_attn_heads)
        self.d_ff = getattr(args, 'attn_enc_filter_size', args.encoder_filter_size)
        self.attention_dropout = getattr(args, 'attn_dropout', args.dropout)
        self.rel_pos_buckets = args.rel_pos_buckets

        encoder_positional_encoding = getattr(args, 'encoder_positional_encoding', 'none')
        if encoder_positional_encoding == "transformer":
            self.encoder_pe = PositionalEncoding(
                dropout=args.dropout,
                dim=self.d_model,
                max_len=1024
            )
        else:
            self.encoder_pe = None
            self._dropout = nn.Dropout(p=args.dropout)

        rel_pos = getattr(args, 'rel_pos', 'emb_only')
        if rel_pos in ["enc_only", "emb_only"]:
            self.u = nn.Parameter(torch.randn(self.d_model), requires_grad=True)
            self.v = nn.Parameter(torch.randn(self.d_model), requires_grad=True)
        else:
            self.u = None
            self.v = None

        shared_encoder_layer = getattr(args, 'shared_encoder_layer', 0)
        if shared_encoder_layer == 1:
            self.attention_layer = SALayerXL(
                args, self.d_model, self.heads, self.d_ff, args.dropout, self.attention_dropout,
                self.rel_pos_buckets, self.u, self.v)
        else:
            self.attention_layers = nn.ModuleList(
                [SALayerXL(
                    args, self.d_model, self.heads, self.d_ff, args.dropout, self.attention_dropout,
                    self.rel_pos_buckets, self.u, self.v)
                 for i in range(self.num_layers)])
        self.layer_norm = nn.LayerNorm(self.d_model, eps=1e-6)

    def forward(self, src, lengths, distances):
        if self.encoder_pe is not None:
            emb = self.encoder_pe(src)
            out = emb.transpose(0, 1).contiguous()
        else:
            out = src.transpose(0, 1).contiguous()
            encoder_emb_scale = getattr(self.args, 'encoder_emb_scale', None)
            if encoder_emb_scale == "sqrt":
                out = out * math.sqrt(self.d_model)
            out = self._dropout(out)

        mask = ~sequence_mask(lengths).unsqueeze(1)

        shared_encoder_layer = getattr(self.args, 'shared_encoder_layer', 0)
        if shared_encoder_layer == 1:
            layer = self.attention_layer
            for i in range(self.num_layers):
                out = layer(out, mask, distances)
        else:
            for layer in self.attention_layers:
                out = layer(out, mask, distances)
        out = self.layer_norm(out)

        return out.transpose(0, 1).contiguous()
