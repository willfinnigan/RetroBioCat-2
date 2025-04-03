from typing import Sequence, List

import pandas as pd

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_bkms import does_bkms_exist, download_bkms_model
from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentDataQuery
from rbc2.precedent_identification.similarity_tools import get_fingerprints

data_folder = f'{path_to_data_folder}/bkms'

class BKMS_DataQuery(PrecedentDataQuery):
    df = None
    fp_df = None

    def __init__(self):
        self.product_column: str = 'products'
        self.substrate_columns: List[str] = []
        self.enzyme_column: str = 'EC_Number'
        self.id_column: str = 'ID'

        if does_bkms_exist() == False:
            download_bkms_model()


    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_hdf(f"{data_folder}/bkms_metadata.hdf")

    def query_data(self, query_ids: Sequence[int]=()) -> pd.DataFrame:
        if len(query_ids) == 0:
            raise Exception('Must specify query_ids to get bkms data')

        self.load_df()
        p_df = self.df.loc[query_ids]
        p_df['rxn_smi'] = p_df['reaction_smiles']

        return p_df

    def get_fps(self, smis):
        return get_fingerprints(smis)

if __name__ == '__main__':
    bkms_data = BKMS_DataQuery()
    bkms_data.load_df()
    query_ids = [4417, 4418, 4419, 4420, 4421, 4422, 4423, 4424, 4425, 4426, 4427, 4428, 4429, 4430, 4431, 4432, 4433, 4434, 4435, 4436, 4437, 4438, 4439, 4440, 4441, 4442, 4443, 4444, 4445, 4446, 4447, 4448, 4449, 4450, 4451, 4452, 4453, 4454, 4455, 4456, 4457, 4458]
    result = bkms_data.query_data(query_ids)
    print(result.iloc[0])
    #print(dict(result.iloc[0]))


