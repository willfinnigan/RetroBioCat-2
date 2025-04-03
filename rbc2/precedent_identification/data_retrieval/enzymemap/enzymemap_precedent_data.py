from typing import Optional

import pandas as pd

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_enzymemap import does_enzymemap_exist, download_enzymemap
from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentDataQuery, get_rxn_smi_vectorized
from rbc2.precedent_identification.similarity_tools import get_fingerprints

data_folder = f'{path_to_data_folder}/enzymemap/brenda'

TEMPLATE_ID_COL_NAME = 'corrected_template_default_id'
METADATA_HDF = 'enzymemap_brenda_metadata.hdf'

class EnzymeMap_DataQuery(PrecedentDataQuery):
    df = None
    fp_df = None

    def __init__(self):
        self.product_column = 'prod_smiles'
        self.substrate_columns = ['reac_smiles']
        self.enzyme_column = 'ec_num'
        self.id_column = None

        if does_enzymemap_exist() == False:
            download_enzymemap()

    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_hdf(f"{data_folder}/{METADATA_HDF}", key='df')

    def query_data(self, template_id: Optional[int] = None) -> pd.DataFrame:
        if template_id is None:
            raise Exception('Must specify a template id to get enzymemap data')

        self.load_df()
        template_df = self.df[self.df[TEMPLATE_ID_COL_NAME] == str(template_id)].copy()
        #template_df['rxn_smi'] = list(template_df.apply(self.get_rxn_smi, axis=1))
        template_df.loc[:, 'rxn_smi'] = get_rxn_smi_vectorized(template_df,
                                                               product_column=self.product_column,
                                                               substrate_columns=self.substrate_columns)

        return template_df

    def get_fps(self, smis):
        return get_fingerprints(smis)

if __name__ == '__main__':
    em_data = EnzymeMap_DataQuery()
    em_data.load_df()
    result = em_data.query_data(1)
    #print(result)
    #print(dict(result.iloc[0]))
    #print(em_data.df.info())

    #df = pd.read_hdf(f"{data_folder}/{METADATA_HDF}", key='df')

    # convert the corrected_template_default_id column from int64 to a normal int
    #df['corrected_template_default_id'] = df['corrected_template_default_id'].astype(str)

    # save the df as hdf again
    #df.to_hdf(f"{data_folder}/{METADATA_HDF}", key='df', mode='w')


