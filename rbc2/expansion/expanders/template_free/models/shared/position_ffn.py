"""Position feed-forward network from "Attention is All You Need"."""
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.utils import skip_init


class ActivationFunction(object):
    relu = "relu"
    gelu = "gelu"
    silu = "silu"
    gated_gelu = "gated-gelu"


ACTIVATION_FUNCTIONS = {
    ActivationFunction.relu: F.relu,
    ActivationFunction.gelu: F.gelu,
    ActivationFunction.silu: F.silu,
    ActivationFunction.gated_gelu: F.gelu,
}


class PositionwiseFeedForward(nn.Module):
    """A two-layer Feed-Forward-Network with residual layer norm.

    Args:
        d_model (int): the size of input for the first-layer of the FFN.
        d_ff (int): the hidden layer size of the second-layer of the FNN.
        dropout (float): dropout probability in [0, 1).
    """

    def __init__(self, d_model, d_ff, dropout=0.1,
                 activation_fn=ActivationFunction.relu,
                 add_ffnbias=True, parallel_residual=False,
                 norm_eps=1e-6, use_ckpting=[], parallel_gpu=1):
        super(PositionwiseFeedForward, self).__init__()
        self.w_1 = skip_init(
            nn.Linear, in_features=d_model,
            out_features=d_ff // parallel_gpu, bias=add_ffnbias)
        self.w_2 = skip_init(
            nn.Linear, in_features=d_ff // parallel_gpu,
            out_features=d_model, bias=add_ffnbias)
        self.layer_norm = nn.LayerNorm(d_model, eps=norm_eps)
        self.parallel_residual = parallel_residual
        self.dropout_p = dropout
        self.dropout_1 = nn.Dropout(dropout)
        self.dropout_2 = nn.Dropout(dropout)
        self.activation = ACTIVATION_FUNCTIONS[activation_fn]
        if activation_fn == "silu" or activation_fn == "gated-gelu":
            self.w_3 = skip_init(
                nn.Linear, in_features=d_model,
                out_features=d_ff // parallel_gpu, bias=add_ffnbias)
        else:
            self.w_3 = None
        self.parallel_gpu = parallel_gpu

    def forward(self, x):
        """Layer definition.

        Args:
            x: ``(batch_size, input_len, model_dim)``

        Returns:
            (FloatTensor): Output ``(batch_size, input_len, model_dim)``.
        """
        if not self.parallel_residual:
            norm_x = self.layer_norm(x)
        else:
            norm_x = x.clone()
        inter = self.w_1(norm_x)
        inter = self.activation(inter)
        if self.w_3 is not None:
            inter.mul_(self.w_3(norm_x))
        if self.dropout_p > 0:
            inter = self.dropout_1(inter)
        inter = self.w_2(inter)
        if self.dropout_p > 0:
            inter = self.dropout_2(inter)
        return inter + x

    def update_dropout(self, dropout):
        self.dropout_1.p = dropout
        self.dropout_2.p = dropout
