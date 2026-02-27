import torch
import torch.nn as nn
import warnings
from rbc2.expansion.expanders.template_free.models.shared.positional_encoding import PositionalEncoding


class Elementwise(nn.ModuleList):
    """A simple network container that applies modules elementwise."""

    def __init__(self, merge=None, *args):
        assert merge in [None, "first", "concat", "sum", "mlp"]
        self.merge = merge
        super(Elementwise, self).__init__(*args)

    def forward(self, emb):
        emb_ = [feat.squeeze(2) for feat in emb.split(1, dim=2)]
        assert len(self) == len(emb_)
        emb_out = []
        for f, x in zip(self, emb_):
            emb_out.append(f(x))
        if self.merge == "first":
            return emb_out[0]
        elif self.merge == "concat" or self.merge == "mlp":
            return torch.cat(emb_out, 2)
        elif self.merge == "sum":
            return sum(emb_out)
        else:
            return emb_out


class Embeddings(nn.Module):
    """Words embeddings for encoder/decoder.

    Args:
        word_vec_size (int): size of the dictionary of embeddings.
        word_padding_idx (int): padding index for words in the embeddings.
        word_vocab_size (int): size of dictionary of embeddings for words.
        position_encoding (bool): see PositionalEncoding
        dropout (float): dropout probability.
    """

    def __init__(self, word_vec_size, word_vocab_size, word_padding_idx,
                 position_encoding=False, feat_merge="concat",
                 feat_vec_exponent=0.7, feat_vec_size=-1,
                 feat_padding_idx=[], feat_vocab_sizes=[],
                 dropout=0, sparse=False, fix_word_vecs=False):
        self._validate_args(feat_merge, feat_vocab_sizes, feat_vec_exponent,
                            feat_vec_size, feat_padding_idx)

        if feat_padding_idx is None:
            feat_padding_idx = []
        self.word_padding_idx = word_padding_idx
        self.word_vec_size = word_vec_size

        vocab_sizes = [word_vocab_size]
        emb_dims = [word_vec_size]
        pad_indices = [word_padding_idx]

        if feat_merge == 'sum':
            feat_dims = [word_vec_size] * len(feat_vocab_sizes)
        elif feat_vec_size > 0:
            feat_dims = [feat_vec_size] * len(feat_vocab_sizes)
        else:
            feat_dims = [int(vocab ** feat_vec_exponent) for vocab in feat_vocab_sizes]
        vocab_sizes.extend(feat_vocab_sizes)
        emb_dims.extend(feat_dims)
        pad_indices.extend(feat_padding_idx)

        emb_params = zip(vocab_sizes, emb_dims, pad_indices)
        embeddings = [nn.Embedding(vocab, dim, padding_idx=pad, sparse=sparse)
                      for vocab, dim, pad in emb_params]
        emb_luts = Elementwise(feat_merge, embeddings)

        self.embedding_size = (sum(emb_dims) if feat_merge == 'concat'
                               else word_vec_size)

        super(Embeddings, self).__init__()
        self.make_embedding = nn.Sequential()
        self.make_embedding.add_module('emb_luts', emb_luts)

        if feat_merge == 'mlp' and len(feat_vocab_sizes) > 0:
            in_dim = sum(emb_dims)
            mlp = nn.Sequential(nn.Linear(in_dim, word_vec_size), nn.ReLU())
            self.make_embedding.add_module('mlp', mlp)

        self.position_encoding = position_encoding

        if self.position_encoding:
            pe = PositionalEncoding(dropout, self.embedding_size)
            self.make_embedding.add_module('pe', pe)

        if fix_word_vecs:
            self.word_lut.weight.requires_grad = False

    def _validate_args(self, feat_merge, feat_vocab_sizes, feat_vec_exponent,
                       feat_vec_size, feat_padding_idx):
        if feat_merge == "sum":
            if feat_vec_exponent != 0.7:
                warnings.warn("Merging with sum, but got non-default feat_vec_exponent.")
            if feat_vec_size != -1:
                warnings.warn("Merging with sum, but got non-default feat_vec_size.")
        elif feat_vec_size > 0:
            if feat_vec_exponent != -1:
                warnings.warn("Not merging with sum and positive feat_vec_size, "
                              "but got non-default feat_vec_exponent.")
        else:
            if feat_vec_exponent <= 0:
                raise ValueError("Using feat_vec_exponent to determine feature vec size, "
                                 "but got feat_vec_exponent <= 0.")
        n_feats = len(feat_vocab_sizes)
        if n_feats != len(feat_padding_idx):
            raise ValueError("Got unequal number of feat_vocab_sizes and "
                             "feat_padding_idx ({:d} != {:d})".format(
                                n_feats, len(feat_padding_idx)))

    @property
    def word_lut(self):
        return self.make_embedding[0][0]

    @property
    def emb_luts(self):
        return self.make_embedding[0]

    def forward(self, source, step=None):
        """Computes the embeddings for words and features.

        Args:
            source (LongTensor): index tensor ``(len, batch, nfeat)``
        Returns:
            FloatTensor: Word embeddings ``(len, batch, embedding_size)``
        """
        if self.position_encoding:
            for i, module in enumerate(self.make_embedding._modules.values()):
                if i == len(self.make_embedding._modules.values()) - 1:
                    source = module(source, step=step)
                else:
                    source = module(source)
        else:
            source = self.make_embedding(source)
        return source

    def update_dropout(self, dropout):
        if self.position_encoding:
            self._modules['make_embedding'][1].dropout.p = dropout
