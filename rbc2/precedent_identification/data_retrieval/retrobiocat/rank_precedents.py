import pandas as pd

"""
Similar precedents are ranked according to activity.
For each unique similar product which is identified, we get the topn enzymes.
To do this, we first score the catagorical data.  Next rank within each caterogy based on any quantitative data 
(similarity or conversion). Specific activity or conversion data can be set as preferred here.
Once the data is ranked, we take the topn enzymes for each unique similar product.  

"""

def score_categorial(cat):
    if cat == 'high':
        return 4
    elif cat == 'medium':
        return 3
    elif cat == 'low':
        return 2
    else:
        return 0

def score_datatype(row, prefer_spec_act=True):
    if prefer_spec_act == True:
        if pd.isna(row['specific_activity']) == False:
            return 2
        elif pd.isna(row['conversion']) == False:
            return 1
    else:
        if pd.isna(row['conversion']) == False:
            return 2
        elif pd.isna(row['specific_activity']) == False:
            return 1
    return 0

def add_data_scores(df, prefer_spec_act):
    df['score_cat'] = df['categorical'].apply(score_categorial)
    df['score_datatype'] = df.apply(score_datatype, axis=1, prefer_spec_act=prefer_spec_act)
    return df

def best_examples_per_product(df, top_n, prefer_spec_act=True):
    # order by high, medium, low, none
    # order within each category by spec_act > conv if prefer_spec_act=True, otherwise opposite
    # order within those categories by number
    # take top_n

    products = df['product_1_smiles'].unique()
    product_dfs = []
    for smi in products:
        p_df = df[df['product_1_smiles'] == smi]

        # sort data depending on whether or not we prefer specific activity over conversion
        if prefer_spec_act == True:
            p_df = p_df.sort_values(['score_cat', 'score_datatype', 'specific_activity', 'conversion'],
                             ascending=[False, False, False, False])
        else:
            p_df = p_df.sort_values(['score_cat', 'score_datatype', 'conversion', 'specific_activity'],
                             ascending=[False, False, False, False])

        # drop duplicate enzyme names so we only get 1 enzyme entry (encase there are multiple).
        # The top (best) will be kept.
        p_df = p_df.drop_duplicates(subset=['enzyme_name'], keep='first')

        top_df = p_df[0:top_n]
        product_dfs.append(top_df)

    return pd.concat(product_dfs)

def get_best_enzymes(df: pd.DataFrame, top_n: int, prefer_spec_act=True) -> pd.DataFrame:
    if df.shape[0] == 0:
        return df

    add_data_scores(df, prefer_spec_act=prefer_spec_act)
    enzyme_types = df['enzyme_type'].unique()
    group_dfs = []
    for et in enzyme_types:
        et_df = df[df['enzyme_type'] == et]
        best_et_df = best_examples_per_product(et_df, top_n, prefer_spec_act=prefer_spec_act)
        group_dfs.append(best_et_df)

    return pd.concat(group_dfs)

if __name__ == '__main__':
    test_data = {'product_1_smiles': ['CCC=O', 'CCC=O', 'CCC=O', 'CCC=O', 'CCC=O', 'CCCCN'],
                 'categorical': ['High', 'Medium', 'Medium', 'Medium', None, 'High'],
                 'enzyme_type': ['CAR', 'CAR', 'CAR', 'CAR', 'CAR', 'CAR'],
                 'enzyme_name': ['mpCAR', 'niCAR', 'srCAR', 'noCAR', 'tpCAR', 'mpCAR'],
                 'conversion': [None, None, 50, 45, None, None],
                 'specific_activity': [1.0, 0.7, None, None, None, 1.0]}

    df = pd.DataFrame(test_data)

    df = get_best_enzymes(df, 2, prefer_spec_act=False)
    for i, row in df.iterrows():
        print(dict(row))
        print()



