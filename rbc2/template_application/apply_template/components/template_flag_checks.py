

# True then filter
from typing import List

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def is_reaction_cyclic(target_smi: str, substrate_smis: List[str]) -> bool:
    """ If any of the reaction substrates matches the reaction reactant,
    then the reaction is cyclic"""
    for smi in substrate_smis:
        if smi == target_smi:
            return True
    return False

def does_reaction_keep_same_number_of_rings(target_smi: str, substrate_smis: List[str]) -> bool:
    mol_r = Chem.MolFromSmiles(substrate_smis[0])
    mol_p = Chem.MolFromSmiles(target_smi)
    change = rdMolDescriptors.CalcNumRings(mol_p) - rdMolDescriptors.CalcNumRings(mol_r)
    if change >= 1:
        return False
    return True

def is_intramolecular(substrate_smis: List[str]) -> bool:
    return len(substrate_smis) == 1

def is_dimer(substrate_smis: List[str]) -> bool:
    if (len(set(substrate_smis)) != 1):  # or len(substrate_smis) != 2  - current retrobiocat code gets sets
        return False  # len of set is not 1, or num substrates is not 2, then it's *not* be a dimer
    return True
