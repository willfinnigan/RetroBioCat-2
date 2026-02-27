"""Data utilities for template-free models. Adapted from ASKCOS for inference only."""

import logging
import re
import numpy as np
import torch
from rdkit import Chem

from rbc2.expansion.expanders.template_free.utils.chem_utils import (
    get_atom_features_sparse, get_bond_features_sparse,
    ATOM_FDIMS, BOND_FDIMS
)
from rbc2.expansion.expanders.template_free.utils.distance_calculator import DistanceCalculator

logger = logging.getLogger(__name__)

_distance_calculator = DistanceCalculator()


def load_vocab_from_file(vocab_file: str) -> dict:
    """Load vocabulary from a text file."""
    vocab = {}
    with open(vocab_file, "r") as f:
        for i, line in enumerate(f):
            token = line.strip().split("\t")[0]
            vocab[token] = i
    return vocab


def tokenize_smi(smi: str) -> str:
    """Tokenize a SMILES string."""
    pattern = r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == "".join(tokens), f"Tokenization failed for {smi}"
    return " ".join(tokens)


def get_egraph_features_from_smi(smi: str) -> tuple:
    """Get expanded graph features from a SMILES string for inference.

    Returns:
        Tuple of (node_features, node_features_len, adjacencies, adjacencies_len)
    """
    node_features = []
    adjacencies = []

    smi = "".join(smi.split())
    if not smi.strip():
        smi = "CC"

    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumBonds() == 0:
        mol = Chem.MolFromSmiles("CC")

    for a in mol.GetAtoms():
        node_feature = [1]  # type = atom
        atom_feature = get_atom_features_sparse(a)
        offset = 3  # 0 reserved, 1 for atom, 2 for bond
        for i in range(len(ATOM_FDIMS)):
            atom_feature[i] += offset
            offset += ATOM_FDIMS[i]
        node_feature.extend(atom_feature)
        node_feature.extend([0, 0, 0, 0])
        node_features.append(node_feature)
    n_atoms = len(node_features)

    for bond_index, b in enumerate(mol.GetBonds()):
        node_feature = [2]  # type = bond
        node_feature.extend([0] * len(ATOM_FDIMS))
        offset = 3 + sum(ATOM_FDIMS)
        bond_feature = get_bond_features_sparse(b)
        for i in range(len(BOND_FDIMS)):
            bond_feature[i] += offset
            offset += BOND_FDIMS[i]
        node_feature.extend(bond_feature)
        node_features.append(node_feature)

        adjacencies.append([bond_index + n_atoms, b.GetBeginAtomIdx()])
        adjacencies.append([bond_index + n_atoms, b.GetEndAtomIdx()])

    node_features = np.array(node_features, dtype=np.int32)
    node_features[node_features >= 9999] = 0
    adjacencies = np.array(adjacencies, dtype=np.int32)

    return node_features, adjacencies


class EG2RSBatch:
    """Batch container for EGraph2RSMILES inference."""

    def __init__(self, fnode, node_lengths, distances,
                 tgt_token_ids=None, tgt_lengths=None):
        self.fnode = fnode
        self.node_lengths = node_lengths
        self.distances = distances
        self.tgt_token_ids = tgt_token_ids
        self.tgt_lengths = tgt_lengths
        self.size = len(node_lengths)

    def to(self, device):
        self.fnode = self.fnode.to(device)
        self.node_lengths = self.node_lengths.to(device)
        if self.distances is not None:
            self.distances = self.distances.to(device)
        if self.tgt_token_ids is not None:
            self.tgt_token_ids = self.tgt_token_ids.to(device)
        if self.tgt_lengths is not None:
            self.tgt_lengths = self.tgt_lengths.to(device)


def create_batch_from_smi(smi: str, rel_pos_buckets: int = 31,
                          device: torch.device = None) -> EG2RSBatch:
    """Create an EG2RSBatch from a single SMILES string for inference.

    Args:
        smi: Product SMILES string.
        rel_pos_buckets: Number of relative position buckets.
        device: torch device to place tensors on.

    Returns:
        EG2RSBatch ready for model input.
    """
    node_features, adjacencies = get_egraph_features_from_smi(smi)
    node_features_len = node_features.shape[0]
    max_len = node_features_len

    # Build adjacency matrix and compute distances
    adjacency_matrix = np.zeros((node_features_len, node_features_len), dtype=np.bool_)
    for a, b in adjacencies:
        adjacency_matrix[a, b] = True
        adjacency_matrix[b, a] = True

    distance = _distance_calculator.calculate(
        adjacency_matrix, node_features_len, max_distance=node_features_len)

    # Bucket distances
    distance[distance >= (rel_pos_buckets - 1)] = rel_pos_buckets - 1
    distance[distance == 0] = rel_pos_buckets
    np.fill_diagonal(distance, 0)

    # Convert to tensors (batch size = 1)
    fnode = torch.as_tensor(node_features, dtype=torch.long).unsqueeze(0)
    node_lengths = torch.tensor([node_features_len], dtype=torch.long)
    distances = torch.as_tensor(distance, dtype=torch.long).unsqueeze(0)

    batch = EG2RSBatch(
        fnode=fnode,
        node_lengths=node_lengths,
        distances=distances)

    if device is not None:
        batch.to(device)

    return batch
