from typing import List

import pandas as pd

from rbc2.configs.data_path import path_to_data_folder
from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentDataQuery, get_rxn_smi_vectorized
from rbc2.precedent_identification.similarity_tools import make_fp_df, get_fingerprints_from_fpdf

data_folder = f'{path_to_data_folder}/retrobiocat'

class RetroBioCatLocalPrecedentDataQuery(PrecedentDataQuery):
    df = None
    fp_df = None

    def __init__(self):
        self.product_column: str = 'product_1_smiles'
        self.substrate_columns: List[str] = ['substrate_1_smiles', 'substrate_2_smiles']
        self.enzyme_column: str = 'enzyme_type'
        self.id_column: str = '_id'

    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_excel(f"{data_folder}/trial_activity_data.xlsx")

    @classmethod
    def load_fp_df(cls):
        if cls.fp_df is None:
            cls.load_df()
            cls.fp_df = make_fp_df(cls.df, 'product_1_smiles')

    def add_rxn_smis(self):
        if 'rxn_smi' not in self.df.columns:
            # self.df.loc[:, 'rxn_smi'] = self.df.apply(self.get_rxn_smi, axis=1)
            self.df.loc[:, 'rxn_smi'] = get_rxn_smi_vectorized(self.df,
                                                               product_column=self.product_column,
                                                               substrate_columns=self.substrate_columns)



    def get_fp_df(self):
        self.load_fp_df()
        return self.fp_df

    def get_single_fp(self, smi: str):
        self.load_df()
        self.load_fp_df()
        try:
            return self.fp_df.loc[smi]['fp']
        except:
            return None

    def get_fps(self, smis: List[str]):
        self.load_df()
        self.load_fp_df()

        fps, converted_smis = get_fingerprints_from_fpdf(smis, self.fp_df)
        return fps, converted_smis


    def query_data(self, reaction_name: str = 'All', enzyme_types: List[str] = ()) -> pd.DataFrame:
        self.load_df()
        self.add_rxn_smis()
        df = self.df

        if reaction_name != 'All':
            df = df[df['reaction'] == reaction_name]

        if len(enzyme_types) != 0 and 'All' not in enzyme_types:
            df = df[df['enzyme_type'].isin(enzyme_types)]

        return df




if __name__ == '__main__':
    local_data = RetroBioCatLocalPrecedentDataQuery()
    #df = local_data.query_data([], [])
    #print(df.head())
    fps, smis = local_data.get_fps(['CCC=O', 'CCCCC=O'])




    # df = local_data.query_data(['Carboxylic acid reduction'], [])
    # print(df.head())
    #
    # df = local_data.query_data(['Reductive amination'], [])
    # print(df.head())
    #
    # df = local_data.query_data(['All'], ['CAR'])
    # print(df.head())
    #
    # from rbc2.mongo.default_connection import make_default_connection
    # make_default_connection()
    # from rbc2.mongo.model_queries.specificity_data_query import query_specificity_data
    #
    # result = query_specificity_data(['Carboxylic acid reduction'], [], only_reviewed=True)
    # print(result)