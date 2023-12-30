import pytest

from rbc2.reaction_evaluation.starting_material_evaluator import StartingMaterialEvaluator


def test_can_lookup_molecule():
    sme = StartingMaterialEvaluator()
    smi = 'CCC=O'
    available, info = sme.eval(smi)
    assert available == 1

def test_random_molecule_not_available():
    sme = StartingMaterialEvaluator()
    smi = 'not_a_molecule'
    available, info = sme.eval(smi)
    assert available == 0

def test_info_dict_is_returned_with_available():
    sme = StartingMaterialEvaluator()
    smi = 'CCC=O'
    available, info = sme.eval(smi)
    assert 'alfa' in info

def test_can_query_metabolism():
    sme = StartingMaterialEvaluator()
    smi = 'CCC(=O)C(=O)O'
    sme.config.source_mol_mode = 'metabolites'
    available, info = sme.eval(smi)
    assert available == 1
    assert 'alfa' not in info

def test_can_add_custom_smis_to_starting_materials():
    sme = StartingMaterialEvaluator()

    available, info = sme.eval('not_a_molecule')
    assert available == 0

    sme.custom_smiles = ['not_a_molecule']

    available, info = sme.eval('not_a_molecule')
    assert available == 1


@pytest.mark.parametrize('allowed,expected', [[True, 1], [False, 0]])
def test_chiral_molecules_are_found_or_not_when_banned(allowed, expected):
    sme = StartingMaterialEvaluator()
    sme.config.source_mols_can_be_chiral = allowed
    smi = 'C[C@H](N)c1ccccc1'
    available, info = sme.eval(smi)
    assert available == expected