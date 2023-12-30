import pandas as pd
from rbc2.precedent_identification.data_retrieval.retrobiocat.rank_precedents import get_best_enzymes

test_data = {'product_1_smiles': ['CCC=O', 'CCC=O', 'CCC=O', 'CCC=O', 'CCC=O', 'CCCCN'],
             'categorical': ['High', 'Medium', 'Medium', 'Medium', None, 'High'],
             'enzyme_type': ['CAR', 'CAR', 'CAR', 'CAR', 'CAR', 'CAR'],
             'enzyme_name': ['mpCAR', 'niCAR', 'srCAR', 'noCAR', 'tpCAR', 'mpCAR'],
             'conversion': [None, None, 50, 45, None, None],
             'specific_activity': [1.0, 0.7, None, None, None, 1.0]}

def test_top_result_has_categorical_high():
    df = pd.DataFrame(test_data)
    df = get_best_enzymes(df, 2)
    assert df.iloc[0]['categorical'] == 'High'

def test_2nd_result_has_categorical_medium_and_has_specific_activity():
    df = pd.DataFrame(test_data)
    df = get_best_enzymes(df, 2)
    assert df.iloc[1]['categorical'] == 'Medium'
    assert df.iloc[1]['specific_activity'] == 0.7

def test_3rd_result_is_for_another_substrate():
    df = pd.DataFrame(test_data)
    df = get_best_enzymes(df, 2)
    assert df.iloc[2]['product_1_smiles'] != df.iloc[0]['product_1_smiles']

def test_no_data_returns_nothing():
    df = pd.DataFrame({})
    df = get_best_enzymes(df, 2)
    assert df.shape[0] == 0