from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np


def get_mol_fingerprint(rd_mol, radius=2, nBits=2048):
    """
    Returns the Morgan fingerprint of the molecule
    """

    bitvect = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius, nBits=nBits)
    #array = np.zeros((1,))
    array = np.empty(nBits, dtype='float32')
    DataStructs.ConvertToNumpyArray(bitvect, array)
    return array


def get_reaction_fingerprint(product_mol, substrate_mols, radius=2, nBits=2048):
    """
    Returns a difference fingerprint
    """

    product_fp = get_mol_fingerprint(product_mol)
    reactants_fp = sum(get_mol_fingerprint(mol, radius=radius, nBits=nBits) for mol in substrate_mols)

    return product_fp - reactants_fp

