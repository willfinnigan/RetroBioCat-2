"""Graph2SMILES molecular featurization and batch construction.

Contains atom/bond feature extraction (one-hot encoding), directed graph
construction, and batch creation for Graph2SMILES inference.

Note: The molecular vocabulary constants (ATOM_LIST, DEGREES, etc.) are shared
with chem_utils.py (used by EGraph2RSMILES), but NUM_Hs_DICT intentionally
differs: G2S uses {0:0, 1:1, 3:2, 4:3, 5:4} to match its training data,
while EGraph2RSMILES uses {0:0, 1:1, 2:2, 3:3, 4:4}.
"""

import numpy as np
import torch
from rdkit import Chem

from rbc2.expansion.expanders.template_free.utils.distance_calculator import DistanceCalculator

# --- Atom and bond vocabulary (Graph2SMILES variant) ---

ATOM_LIST = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe',
             'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti',
             'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb',
             'W', 'Ru', 'Nb', 'Re', 'Te', 'Rh', 'Ta', 'Tc', 'Ba', 'Bi', 'Hf', 'Mo', 'U', 'Sm', 'Os', 'Ir',
             'Ce', 'Gd', 'Ga', 'Cs', '*', 'unk']
ATOM_DICT = {symbol: i for i, symbol in enumerate(ATOM_LIST)}

DEGREES = list(range(10))
HYBRIDIZATION = [Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
                 Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
                 Chem.rdchem.HybridizationType.SP3D2]
HYBRIDIZATION_DICT = {hb: i for i, hb in enumerate(HYBRIDIZATION)}
FC_DICT = {-1: 0, -2: 1, 1: 2, 2: 3, 0: 4}
VALENCE_DICT = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6}
# Intentionally different from chem_utils.py - matches G2S training data
NUM_Hs_DICT = {0: 0, 1: 1, 3: 2, 4: 3, 5: 4}
CHIRAL_TAG = [Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
              Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
              Chem.rdchem.ChiralType.CHI_UNSPECIFIED]
CHIRAL_TAG_DICT = {ct: i for i, ct in enumerate(CHIRAL_TAG)}
RS_TAG_DICT = {"R": 0, "S": 1, "None": 2}

BOND_TYPES = [None, Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE,
              Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]
BOND_STEREO = [Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOZ,
               Chem.rdchem.BondStereo.STEREONONE]

ATOM_FDIM = [len(ATOM_LIST), len(DEGREES), 5, len(HYBRIDIZATION), 7,
             5, len(CHIRAL_TAG), 3, 2]
BOND_FDIM = 9  # 4 bond types + 3 stereo + conjugated + ring


# --- Feature extraction ---

def get_atom_features(atom):
    """Get one-hot encoded atom features for Graph2SMILES."""
    symbol_id = ATOM_DICT.get(atom.GetSymbol(), ATOM_DICT["unk"])
    if atom.GetSymbol() in ["*", "unk"]:
        features = [0] * sum(ATOM_FDIM)
        if symbol_id < len(ATOM_LIST):
            features[symbol_id] = 1
        return features
    degree_id = min(atom.GetDegree(), 9)
    formal_charge_id = FC_DICT.get(atom.GetFormalCharge(), 4)
    hybridization_id = HYBRIDIZATION_DICT.get(atom.GetHybridization(), 4)
    valence_id = VALENCE_DICT.get(atom.GetTotalValence(), 6)
    num_h_id = NUM_Hs_DICT.get(atom.GetTotalNumHs(), 0)
    chiral_tag_id = CHIRAL_TAG_DICT.get(atom.GetChiralTag(), 2)
    rs_tag = atom.GetPropsAsDict().get("_CIPCode", "None")
    rs_tag_id = RS_TAG_DICT.get(rs_tag, 2)
    is_aromatic = int(atom.GetIsAromatic())
    features = []
    for feat_val, feat_dim in zip(
        [symbol_id, degree_id, formal_charge_id, hybridization_id,
         valence_id, num_h_id, chiral_tag_id, rs_tag_id, is_aromatic],
        ATOM_FDIM):
        one_hot = [0] * feat_dim
        if feat_val < feat_dim:
            one_hot[feat_val] = 1
        features.extend(one_hot)
    return features


def get_bond_features(bond):
    """Get one-hot encoded bond features for Graph2SMILES."""
    bt = bond.GetBondType()
    bond_features = [int(bt == bond_type) for bond_type in BOND_TYPES[1:]]
    bs = bond.GetStereo()
    bond_features.extend([int(bs == stereo) for stereo in BOND_STEREO])
    bond_features.extend([int(bond.GetIsConjugated()), int(bond.IsInRing())])
    return bond_features


# --- G2S Batch ---

class G2SBatch:
    """Batch container for Graph2SMILES inference."""

    def __init__(self, fnode, fmess, agraph, bgraph, atom_scope, bond_scope,
                 tgt_token_ids=None, tgt_lengths=None, distances=None):
        self.fnode = fnode
        self.fmess = fmess
        self.agraph = agraph
        self.bgraph = bgraph
        self.atom_scope = atom_scope
        self.bond_scope = bond_scope
        self.tgt_token_ids = tgt_token_ids
        self.tgt_lengths = tgt_lengths
        self.distances = distances
        self.size = 1

    def to(self, device):
        self.fnode = self.fnode.to(device)
        self.fmess = self.fmess.to(device)
        self.agraph = self.agraph.to(device)
        self.bgraph = self.bgraph.to(device)
        if self.distances is not None:
            self.distances = self.distances.to(device)
        if self.tgt_token_ids is not None:
            self.tgt_token_ids = self.tgt_token_ids.to(device)


# --- Graph construction ---

def build_graph_from_mol(mol):
    """Build directed graph features from an RDKit mol object."""
    n_atoms = mol.GetNumAtoms()

    atom_features = []
    for atom in mol.GetAtoms():
        atom_features.append(get_atom_features(atom))

    # Build adjacency dict: for each atom, list of neighbors
    neighbors = [[] for _ in range(n_atoms)]
    for bond in mol.GetBonds():
        u = bond.GetBeginAtom().GetIdx()
        v = bond.GetEndAtom().GetIdx()
        neighbors[u].append(v)
        neighbors[v].append(u)

    edge_dict = {}
    bond_features = []
    a_graphs = [[] for _ in range(n_atoms)]

    for bond in mol.GetBonds():
        u = bond.GetBeginAtom().GetIdx()
        v = bond.GetEndAtom().GetIdx()
        bf = get_bond_features(bond)

        # Forward edge
        eid_fwd = len(bond_features)
        bond_features.append([u, v] + bf)
        edge_dict[(u, v)] = eid_fwd
        a_graphs[v].append(eid_fwd)

        # Backward edge
        eid_bwd = len(bond_features)
        bond_features.append([v, u] + bf)
        edge_dict[(v, u)] = eid_bwd
        a_graphs[u].append(eid_bwd)

    b_graphs = [[] for _ in range(len(bond_features))]
    for i, bf in enumerate(bond_features):
        u, v = bf[:2]
        for w in neighbors[u]:
            if w != v and (w, u) in edge_dict:
                b_graphs[i].append(edge_dict[(w, u)])

    return atom_features, bond_features, a_graphs, b_graphs


def create_g2s_batch_from_smi(smi: str, rel_pos_buckets: int = 31,
                               compute_graph_distance: bool = True,
                               device: torch.device = None) -> G2SBatch:
    """Create a G2SBatch from a single SMILES string for inference."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumBonds() == 0:
        mol = Chem.MolFromSmiles("CC")

    atom_features, bond_features, a_graphs, b_graphs = build_graph_from_mol(mol)

    n_atoms = len(atom_features)
    n_bonds = len(bond_features)
    max_nn = 11

    # Padding atom (index 0) uses '*' atom features
    star_features = [0] * sum(ATOM_FDIM)
    star_features[ATOM_DICT['*']] = 1
    atom_features = [star_features] + atom_features

    # Bond padding entry + offset atom indices by +1
    bond_features = [[0, 0] + [0] * BOND_FDIM] + bond_features
    for i in range(1, len(bond_features)):
        bond_features[i][0] += 1
        bond_features[i][1] += 1

    # Offset real edge indices by +1, pad with 0s (pointing to dummy entry)
    a_graphs_padded = [[0] * max_nn]
    for ag in a_graphs:
        row = [idx + 1 for idx in ag]
        row.extend([0] * (max_nn - len(row)))
        a_graphs_padded.append(row)
    b_graphs_padded = [[0] * max_nn]
    for bg in b_graphs:
        row = [idx + 1 for idx in bg]
        row.extend([0] * (max_nn - len(row)))
        b_graphs_padded.append(row)

    fnode = torch.tensor(atom_features, dtype=torch.float)
    fmess = torch.tensor(bond_features, dtype=torch.float)
    agraph = torch.tensor(a_graphs_padded, dtype=torch.long)
    bgraph = torch.tensor(b_graphs_padded, dtype=torch.long)

    # Compute distances
    distances = None
    if compute_graph_distance and n_atoms > 0:
        adj = Chem.rdmolops.GetAdjacencyMatrix(mol).astype(np.bool_)
        calc = DistanceCalculator()
        dist = calc.calculate(adj, n_atoms, max_distance=n_atoms)
        dist[dist >= (rel_pos_buckets - 1)] = rel_pos_buckets - 1
        dist[dist == 0] = rel_pos_buckets
        np.fill_diagonal(dist, 0)
        distances = torch.tensor(dist, dtype=torch.long).unsqueeze(0)

    # atom_scope starts at 1 (offset for padding entry)
    batch = G2SBatch(
        fnode=fnode, fmess=fmess, agraph=agraph, bgraph=bgraph,
        atom_scope=[[(1, n_atoms)]], bond_scope=[[(1, n_bonds)]],
        distances=distances)

    if device is not None:
        batch.to(device)

    return batch
