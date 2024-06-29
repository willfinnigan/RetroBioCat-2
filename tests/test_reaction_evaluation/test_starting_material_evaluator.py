import pytest
from rbc2 import EcoliSME, CommercialSME


def test_can_lookup_molecule():
    sme = CommercialSME()
    smi = 'CCC=O'
    available, info = sme.eval(smi)
    assert available == True

def test_random_molecule_not_available():
    sme = CommercialSME()
    smi = 'not_a_molecule'
    available, info = sme.eval(smi)
    assert available == False

def test_info_dict_is_returned_with_available():
    sme = CommercialSME()
    smi = 'CCC=O'
    available, info = sme.eval(smi)
    assert 'alfa' in info

def test_can_query_metabolism():
    sme = EcoliSME()
    smi = 'CCC(=O)C(=O)O'
    available, info = sme.eval(smi)
    assert available == True
    assert 'alfa' not in info

def test_can_add_custom_smis_to_starting_materials():
    sme = CommercialSME()
    available, info = sme.eval('not_a_molecule')
    assert available == False

    sme = CommercialSME()
    sme.custom_smiles = ['not_a_molecule']
    available, info = sme.eval('not_a_molecule')
    assert available == True


@pytest.mark.parametrize('allowed,expected', [[True, 1], [False, 0]])
def test_chiral_molecules_are_found_or_not_when_banned(allowed, expected):
    sme = CommercialSME()
    sme.source_mols_can_be_chiral = allowed
    smi = 'C[C@H](N)c1ccccc1'
    available, info = sme.eval(smi)
    assert available == expected