import torch
import torch.nn as nn
from rbc2.expansion.expanders.template_free.models.graph2smiles.dgat import DGATEncoder
from rbc2.expansion.expanders.template_free.models.graph2smiles.dgcn import DGCNEncoder
from typing import Tuple


class GraphFeatEncoder(nn.Module):
    """GraphFeatEncoder encodes molecules by using features of atoms and bonds."""

    def __init__(self, args, n_atom_feat: int, n_bond_feat: int):
        super().__init__()
        self.args = args
        self.n_atom_feat = n_atom_feat
        self.n_bond_feat = n_bond_feat

        if args.mpn_type == "dgcn":
            MPNClass = DGCNEncoder
        elif args.mpn_type == "dgat":
            MPNClass = DGATEncoder
        else:
            raise NotImplementedError(f"Unsupported mpn_type: {args.mpn_type}!")

        self.mpn = MPNClass(
            args,
            input_size=n_atom_feat + n_bond_feat,
            node_fdim=n_atom_feat)

    def forward(self, reaction_batch) -> Tuple[torch.Tensor, None]:
        fnode = reaction_batch.fnode
        fmess = reaction_batch.fmess
        agraph = reaction_batch.agraph
        bgraph = reaction_batch.bgraph

        hnode = fnode.clone()
        fmess1 = hnode.index_select(index=fmess[:, 0].long(), dim=0)
        fmess2 = fmess[:, 2:].clone()
        hmess = torch.cat([fmess1, fmess2], dim=-1)

        hatom, _ = self.mpn(hnode, hmess, agraph, bgraph, mask=None)
        return hatom, None
