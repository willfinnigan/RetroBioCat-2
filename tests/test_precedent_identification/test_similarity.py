import pandas as pd

from rbc2.precedent_identification.similarity_tools import get_single_fp, bulk_similarity, get_fingerprints, \
    get_fingerprints_from_fpdf, make_fp_df


def test_bulk_similarity():
    target_smi = 'CCCC=O'
    compare_with = ['CCCO', 'CCCCN', 'CCC(C)CC(=O)O']
    target_fp = get_single_fp(target_smi)
    sims = bulk_similarity(target_fp, compare_with, get_fingerprints)
    assert sims == {"CCCO": 0.2, 'CCCCN': 0.2727272727272727, 'CCC(C)CC(=O)O': 0.25}

def test_bulk_similarity_empty_string():
    target_smi = 'CCCCCCCCC=O'
    compare_with = ['', '', '']
    target_fp = get_single_fp(target_smi)
    sims = bulk_similarity(target_fp, compare_with, get_fingerprints)
    assert sims == {"": 0}

def test_get_fingerprints():
    smis = ['CCCO', 'CCCCN', 'CCC(C)CC(=O)O', 'nan', None]
    fps, fp_smis = get_fingerprints(smis)
    assert len(fps) == len(fp_smis)

def test_make_fp_df():
    smis = ['CCCO', 'CCCCN', 'CCC(C)CC(=O)O']
    df = pd.DataFrame({'smiles': smis})
    fp_df = make_fp_df(df, 'smiles')
    assert fp_df.loc['CCCO'][0] == get_single_fp('CCCO')
    assert fp_df.shape[0] == 3

def test_can_lookup_fp_from_fpdf():
    smis = ['CCCO', 'CCCCN', 'CCC(C)CC(=O)O']
    df = pd.DataFrame({'smiles': smis})
    fp_df = make_fp_df(df, 'smiles')
    fps, c_smis = get_fingerprints_from_fpdf(['CCCO'], fp_df)
    assert fps == [get_single_fp('CCCO')]
