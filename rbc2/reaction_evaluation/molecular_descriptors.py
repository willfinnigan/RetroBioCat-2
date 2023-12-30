from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt


def get_mw(smi):
    mol = Chem.MolFromSmiles(smi)
    return ExactMolWt(mol)