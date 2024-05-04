from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def num_rings(outcomes, target):
    """Given a reaction SMILES, obtains the change in the number of rings

    Parameters:
        reaction (SMILES): The reaction without agents in the form of SMILES

    Returns:
        num_rings (int): The change in the number of rings in the reaction
    """
    reaction_smi = '>>'.join([outcomes[0], target])

    mol_r = Chem.MolFromSmiles(reaction_smi.split('>>')[0])
    mol_p = Chem.MolFromSmiles(reaction_smi.split('>>')[-1])

    return rdMolDescriptors.CalcNumRings(mol_p) - rdMolDescriptors.CalcNumRings(mol_r)

