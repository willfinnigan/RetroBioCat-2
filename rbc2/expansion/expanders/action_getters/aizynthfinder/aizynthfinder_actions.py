import time

from rdkit import Chem
import numpy as np
import pandas as pd

from rbc2.configs.download_data_files.download_aizynthfinder import does_aizynthfinder_exist, \
    download_aizynthfinder_model
from rbc2.utils.add_logger import add_logger
from rbc2.configs.data_path import path_to_data_folder

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.utils import load_keras_models, fingerprints

data_folder = f'{path_to_data_folder}/aizynthfinder'

class AizynthfinderActionGetter():

    def __init__(self,
                 template_column='retro_template',
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 log_level='WARNING'):

        self.logger = add_logger('AIZynthfinder_Actions', level=log_level)
        self.policy_model = None
        self.templates = None

        self.template_column = template_column
        self.cutoff_cumulative = cutoff_cumulative
        self.cutoff_number = cutoff_number

        if does_aizynthfinder_exist() == False:
            download_aizynthfinder_model()

    def load_model(self):
        if self.policy_model == None:
            policy_path = data_folder + '/uspto_model.hdf5'
            self.policy_model = load_keras_models.LocalKerasModel(policy_path)
        if self.templates == None:
            templates_path = data_folder + '/uspto_templates.hdf5'
            self.templates = pd.read_hdf(templates_path, "table")

    def get_actions(self, smi):
        reactions = []
        priors = []
        template_column = self.template_column

        mol = Chem.MolFromSmiles(smi)

        all_transforms_prop = self._predict(mol)

        probable_transforms_idx = self._cutoff_predictions(all_transforms_prop)

        possible_moves = self.templates.iloc[probable_transforms_idx]
        probs = all_transforms_prop[probable_transforms_idx]

        priors.extend(probs)
        for idx, (move_index, move) in enumerate(possible_moves.iterrows()):
            metadata = dict(move)
            del metadata[template_column]
            metadata["policy_probability"] = round(float(probs[idx]), 5)
            metadata["template_code"] = move_index

            reaction = {'smarts': move[template_column],
                        'metadata': metadata,
                        'prior': priors[idx]}

            reactions.append(reaction)

        return reactions

    def get_rxns(self, smile):
        if self.policy_model == None:
            self.load_model()

        reactions = self.get_actions(smile)
        rxns = {}
        metadata = {}

        for reaction in reactions:
            name = f"Chem_{reaction['metadata']['classification']}"
            num = 1
            extra_string = f"__{num}"
            while name+extra_string in rxns:
                extra_string = f"__{num}"
                num += 1
            name = name+extra_string
            smarts = reaction['smarts']
            if self._does_smarts_only_one_reactants(smarts):
                rxns[name] = [smarts]
            else:
                rxns[name] = []
            metadata[name] = reaction['metadata']
        return rxns, metadata

    def _predict(self, mol):
        fingerprint = fingerprints.get_mol_fingerprint(mol, 2, nBits=len(self.policy_model))
        fp_arr = fingerprint.reshape([1, len(self.policy_model)])
        return np.array(self.policy_model.predict(fp_arr)).flatten()

    @staticmethod
    def _does_smarts_only_one_reactants(smarts):
        if '>>' not in smarts:
            return False
        if '.' in smarts.split('>>')[0]:
            return False
        return True

    def _cutoff_predictions(self, predictions):
        """
        Get the top transformations, by selecting those that have:
            * cumulative probability less than a threshold (cutoff_cumulative)
            * or at most N (cutoff_number)
        """

        sortidx = np.argsort(predictions)[::-1]
        cumsum = np.cumsum(predictions[sortidx])
        if any(cumsum >= self.cutoff_cumulative):
            maxidx = np.argmin(cumsum < self.cutoff_cumulative)
        else:
            maxidx = len(cumsum)

        maxidx = min(maxidx, self.cutoff_number) or 1
        return sortidx[:maxidx]

if __name__ == '__main__':
    smi = 'CC(C)(C)C1=CC=C(C=C1)C(=O)O'
    getter = AizynthfinderActionGetter()
    rxns, metadata = getter.get_rxns(smi)

    times = []
    for i in range(100):
        t0 = time.time()
        rxns, metadata = getter.get_rxns(smi)
        print(rxns)
        t1 = time.time()
        speed = (t1-t0)*1000
        times.append(speed)
        print(speed)

    print(f"Average time: {np.mean(times)} ms")

