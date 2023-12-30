from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tensorflow import keras
from scipy.special import softmax

from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.configs.data_path import path_to_data_folder

try:
    from tensorflow.python.ops.numpy_ops import np_config
    np_config.enable_numpy_behavior()
except:
    pass


data_folder = f'{path_to_data_folder}/bkms'

class BKMS_Action_Getter():


    def __init__(self,
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 allow_multi_product_templates=False,
                 log_level='WARNING'):
        self.logger = add_logger('BKMSActionGetter', level=log_level)
        self.policy_model = None
        self.templates = None

        self.fp_length = 2048
        self.fp_radius = 2

        self.cutoff_cumulative = cutoff_cumulative
        self.cutoff_number = cutoff_number
        self.allow_multi_product_templates = allow_multi_product_templates


    def load_model(self):
        if self.policy_model is None:
            policy_path = data_folder + '/policy_model'
            self.policy_model = keras.models.load_model(policy_path)
        if self.templates is None:
            templates_path = data_folder + '/bkms_templates_only.hdf'
            self.templates: pd.DataFrame = pd.read_hdf(templates_path, key='df')

    def get_rxns(self, smiles):
        self.load_model()

        fp = self._smiles_to_fp(smiles).reshape(1, -1)
        indices, scores = self._predictions_from_model(fp)
        indices, scores = self._filter_by_cumulative_proba(indices, scores)
        indices, scores = self._filter_by_max_number(indices, scores)
        possible_moves = self.templates.iloc[indices]

        smarts_dict, metadata_dict = {}, {}
        for i, (index, row) in enumerate(possible_moves.iterrows()):
            name = f"BKMS_{i+1}"
            smarts = row['reaction_smarts']
            if self._does_template_have_multiple_products(smarts):
                continue
            metadata = dict(row)
            metadata['score'] = float(scores[i])
            smarts_dict[name] = [smarts]
            metadata_dict[name] = metadata

        self.logger.debug(f'{len(smarts_dict)} options retrieved for {smiles}')
        return smarts_dict, metadata_dict

    def _does_template_have_multiple_products(self, smarts):
        if self.allow_multi_product_templates:
            return False

        products = smarts.split('>>')[0].split('.')
        if len(products) != 1:
            return True

        return False

    def _predictions_from_model(self, fp):
        scores = self.policy_model(fp).reshape(-1)
        scores = softmax(scores)
        indices = np.argsort(-scores)
        scores = scores[indices]
        return indices, scores

    def _filter_by_cumulative_proba(self, indices, scores):
        cum_scores = np.cumsum(scores)
        scores = scores[cum_scores <= self.cutoff_cumulative]
        indices = indices[cum_scores <= self.cutoff_cumulative]
        return indices, scores

    def _filter_by_max_number(self, indices, scores):
        indices = indices[:self.cutoff_number]
        scores = scores[:self.cutoff_number]
        return indices, scores

    def _smiles_to_fp(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return np.zeros((self.fp_length,), dtype=np.float32)
        return np.array(
            AllChem.GetMorganFingerprintAsBitVect(mol, self.fp_radius, nBits=self.fp_length, useChirality=True), dtype=np.float32
        )


if __name__ == '__main__':
    getter = BKMS_Action_Getter()
    smarts_dict, metadata_dict = getter.get_rxns('CCC=O')
    print(smarts_dict)
    print(metadata_dict)