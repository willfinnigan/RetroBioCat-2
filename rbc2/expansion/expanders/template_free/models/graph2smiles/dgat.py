import torch
import torch.nn as nn
from rbc2.expansion.expanders.template_free.models.graph2smiles.model_utils import index_select_ND
from typing import Tuple


class DGATGRU(nn.Module):
    """GRU Message Passing layer with attention."""
    def __init__(self, args, input_size: int, h_size: int, depth: int):
        super().__init__()
        self.args = args
        self.input_size = input_size
        self.h_size = h_size
        self.depth = depth
        self._build_layer_components()
        self._build_attention()

    def _build_layer_components(self):
        self.W_z = nn.Linear(self.input_size + self.h_size, self.h_size)
        self.W_r = nn.Linear(self.input_size, self.h_size, bias=False)
        self.U_r = nn.Linear(self.h_size, self.h_size)
        self.W_h = nn.Linear(self.input_size + self.h_size, self.h_size)

    def _build_attention(self):
        self.leaky_relu = nn.LeakyReLU()
        self.head_count = self.args.dgat_attn_heads
        self.dim_per_head = self.h_size // self.head_count
        self.attn_alpha = nn.Parameter(
            torch.Tensor(1, 1, self.head_count, 2 * self.dim_per_head), requires_grad=True)
        self.attn_bias = nn.Parameter(
            torch.Tensor(1, 1, self.head_count), requires_grad=True)
        self.attn_W_q = nn.Linear(self.input_size, self.h_size, bias=True)
        self.attn_W_k = nn.Linear(self.h_size, self.h_size, bias=True)
        self.attn_W_v = nn.Linear(self.h_size, self.h_size, bias=True)
        self.softmax = nn.Softmax(dim=1)
        self.dropout = nn.Dropout(self.args.dropout)
        self.attn_dropout = nn.Dropout(self.args.attn_dropout)

    def GRU(self, x, h_nei):
        n_node, max_nn, h_size = h_nei.size()
        head_count = self.head_count
        dim_per_head = self.dim_per_head

        q = self.attn_W_q(x).unsqueeze(1).repeat(1, max_nn, 1)
        q = q.reshape(n_node, max_nn, head_count, dim_per_head)
        k = self.attn_W_k(h_nei).reshape(n_node, max_nn, head_count, dim_per_head)
        v = self.attn_W_v(h_nei).reshape(n_node, max_nn, head_count, dim_per_head)

        qk = torch.cat([q, k], dim=-1)
        qk = self.leaky_relu(qk)
        attn_score = (qk * self.attn_alpha).sum(dim=-1) + self.attn_bias
        attn_mask = (h_nei.sum(dim=2) == 0).unsqueeze(2)
        attn_score = attn_score.masked_fill(attn_mask, -1e18)
        attn_weight = self.softmax(attn_score).unsqueeze(3)
        attn_context = (attn_weight * v).reshape(n_node, max_nn, h_size)
        sum_h = attn_context.sum(dim=1)

        z_input = torch.cat([x, sum_h], dim=1)
        z = torch.sigmoid(self.W_z(z_input))
        r_1 = self.W_r(x)
        r_2 = self.U_r(sum_h)
        r = torch.sigmoid(r_1 + r_2)
        sum_gated_h = r * sum_h
        h_input = torch.cat([x, sum_gated_h], dim=1)
        pre_h = torch.tanh(self.W_h(h_input))
        new_h = (1.0 - z) * sum_h + z * pre_h
        return new_h

    def forward(self, fmess, bgraph):
        h = torch.zeros(fmess.size()[0], self.h_size, device=fmess.device)
        mask = torch.ones(h.size()[0], 1, device=h.device)
        mask[0, 0] = 0
        for i in range(self.depth):
            h_nei = index_select_ND(h, 0, bgraph)
            h = self.GRU(fmess, h_nei)
            h = h * mask
        return h


class DGATEncoder(nn.Module):
    """MessagePassing Network encoder using DGAT."""
    def __init__(self, args, input_size: int, node_fdim: int):
        super().__init__()
        self.args = args
        self.h_size = args.encoder_hidden_size
        self.depth = args.encoder_num_layers
        self.input_size = input_size
        self.node_fdim = node_fdim
        self.head_count = self.args.dgat_attn_heads
        self.dim_per_head = self.h_size // self.head_count
        self.leaky_relu = nn.LeakyReLU()
        self._build_layers()
        self._build_attention()

    def _build_layers(self):
        self.W_o = nn.Sequential(nn.Linear(self.node_fdim + self.h_size, self.h_size), nn.GELU())
        self.rnn = DGATGRU(self.args, self.input_size, self.h_size, self.depth)

    def _build_attention(self):
        self.attn_alpha = nn.Parameter(
            torch.Tensor(1, 1, self.head_count, 2 * self.dim_per_head), requires_grad=True)
        self.attn_bias = nn.Parameter(
            torch.Tensor(1, 1, self.head_count), requires_grad=True)
        self.attn_W_q = nn.Linear(self.node_fdim, self.h_size, bias=True)
        self.attn_W_k = nn.Linear(self.h_size, self.h_size, bias=True)
        self.attn_W_v = nn.Linear(self.h_size, self.h_size, bias=True)
        self.softmax = nn.Softmax(dim=1)
        self.dropout = nn.Dropout(self.args.dropout)
        self.attn_dropout = nn.Dropout(self.args.attn_dropout)

    def forward(self, fnode, fmess, agraph, bgraph, mask) -> Tuple[torch.Tensor, ...]:
        h = self.rnn(fmess, bgraph)
        nei_message = index_select_ND(h, 0, agraph)

        n_node, max_nn, h_size = nei_message.size()
        head_count = self.head_count
        dim_per_head = self.dim_per_head

        q = self.attn_W_q(fnode).unsqueeze(1).repeat(1, max_nn, 1)
        q = q.reshape(n_node, max_nn, head_count, dim_per_head)
        k = self.attn_W_k(nei_message).reshape(n_node, max_nn, head_count, dim_per_head)
        v = self.attn_W_v(nei_message).reshape(n_node, max_nn, head_count, dim_per_head)

        qk = torch.cat([q, k], dim=-1)
        qk = self.leaky_relu(qk)
        attn_score = (qk * self.attn_alpha).sum(dim=-1) + self.attn_bias
        attn_mask = (nei_message.sum(dim=2) == 0).unsqueeze(2)
        attn_score = attn_score.masked_fill(attn_mask, -1e18)
        attn_weight = self.softmax(attn_score).unsqueeze(3)
        attn_context = (attn_weight * v).reshape(n_node, max_nn, h_size)
        nei_message = attn_context.sum(dim=1)

        node_hiddens = torch.cat([fnode, nei_message], dim=1)
        node_hiddens = self.W_o(node_hiddens)
        if mask is None:
            mask = torch.ones(node_hiddens.size(0), 1, device=fnode.device)
            mask[0, 0] = 0
        return node_hiddens * mask, h
