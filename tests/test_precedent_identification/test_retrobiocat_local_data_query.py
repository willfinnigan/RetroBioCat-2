from rbc2.precedent_identification.data_retrieval.retrobiocat.local_data_query import RetroBioCatLocalPrecedentData
from rbc2.precedent_identification.similarity_tools import get_single_fp, bulk_similarity


def test_local_data_query_object_can_load_data():
    local_data = RetroBioCatLocalPrecedentData()
    local_data.load_df()
    assert local_data.df.shape[0] > 0

def test_local_data_query_object_can_load_fps():
    local_data = RetroBioCatLocalPrecedentData()
    fp_df = local_data.get_fp_df()
    assert fp_df.shape[0] > 0

def test_can_query_local_data():
    local_data = RetroBioCatLocalPrecedentData()

    df = local_data.query_data(reaction_name='Carboxylic acid reduction')
    assert list(df['reaction'].unique()) == ['Carboxylic acid reduction']

    df = local_data.query_data(reaction_name='Reductive amination')
    assert list(df['reaction'].unique()) == ['Reductive amination']

    df = local_data.query_data(enzyme_types=['CAR'])
    assert list(df['enzyme_type'].unique()) == ['CAR']

def test_local_data_query_can_get_fp():
    local_data = RetroBioCatLocalPrecedentData()
    smi = 'OC1CCC(=CBr)CC1'
    fp = local_data.get_single_fp(smi)
    test_fp = get_single_fp(smi)
    assert fp == test_fp

def test_local_data_query_returns_none_if_not_present():
    local_data = RetroBioCatLocalPrecedentData()
    smi = 'C(N)C(N)C(N)C(N)'
    fp = local_data.get_single_fp(smi)
    assert fp == None

def test_local_data_query_can_return_multiple_fps():
    local_data = RetroBioCatLocalPrecedentData()
    smis = ['CCC=O', 'CCCCC=O']
    fps, fp_smis = local_data.get_fps(smis)
    assert len(fps) == len(fp_smis)
    assert len(fps) == len(smis)

def test_sims():
    smis = ['C=CCC(C(=O)NCc1ccccc1)C(C)O', 'CC(C)C[C@@H](C)N', 'CC(O)CCCCO', 'C=CC(O)CCCCC', 'OC1CCOC1',
            'C[C@H](O)c1cccc(Br)c1', 'N[C@@H](C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(Cl)CSC12)c1ccccc1', '',
            'CC1CNC(c2ccccc2)CN1', 'O=C(NCCc1ccccc1)c1ccc2[nH]c3ccccc3c2c1']

    target_smi = 'CCCCCCCCC=O'
    target_fp = get_single_fp(target_smi)

    local_data = RetroBioCatLocalPrecedentData()

    sims = bulk_similarity(target_fp, smis, local_data.get_fps)
    print(sims)


