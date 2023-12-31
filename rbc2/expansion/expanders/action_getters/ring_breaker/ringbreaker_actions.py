import functools
import numpy as np
import pandas as pd

from tensorflow.keras.metrics import top_k_categorical_accuracy
from tensorflow.keras.models import load_model

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataStructs import cDataStructs

from sklearn.preprocessing import LabelEncoder

from rbc2.configs.download_data_files.download_ringbreaker import does_ringbreaker_exist, download_ringbreaker_model
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.configs.data_path import path_to_data_folder

data_folder = f'{path_to_data_folder}/ringbreaker'

class RingBreaker_ActionGetter():


    def __init__(self,
                 cutoff_cumulative=0.995,
                 cutoff_number=10,
                 log_level='WARNING'):

        self.logger = add_logger('RingBreaker_Actions', level=log_level)

        self.cutoff_cumulative = cutoff_cumulative
        self.cutoff_number = cutoff_number

        self.policy_model = None
        self.templates = None
        self.lib = None

        if does_ringbreaker_exist() == False:
            download_ringbreaker_model()

    def init_template_library(self):
        # Initalise template library

        templates_path = f"{data_folder}/data/uspto_ringformations.csv"
        lib = pd.read_csv(templates_path)
        lib = lib.drop(["Unnamed: 0", "index", "selectivity", "outcomes", "ring_change"], axis=1)
        template_labels = LabelEncoder()
        lib['template_code'] = template_labels.fit_transform(lib['template_hash'])
        lib = lib.drop(["reaction_hash", "reactants"], axis=1)

        self.lib = lib.set_index('template_code').T.to_dict('list')

    def init_model(self):
        # Initialise policy
        top10_acc = functools.partial(top_k_categorical_accuracy, k=10)
        top10_acc.__name__ = 'top10_acc'

        top50_acc = functools.partial(top_k_categorical_accuracy, k=50)
        top50_acc.__name__ = 'top50_acc'

        model_path = f"{data_folder}/models/checkpoints/weights.hdf5"
        self.policy = load_model(model_path, custom_objects={'top10_acc': top10_acc, 'top50_acc': top50_acc})

    def load_model(self):
        if self.policy_model is None:
            self.init_model()
        if self.lib is None:
            self.init_template_library()

    def _smiles_to_ecfp(self, product, size=2048):
        mol = Chem.MolFromSmiles(product)
        ecfp = GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((0,), dtype=np.int8)
        cDataStructs.ConvertToNumpyArray(ecfp, arr)
        return arr.reshape([1, 2048])

    def _get_prediction(self, target):
        return self.policy(self._smiles_to_ecfp(target))

    def _get_templates(self, prediction, lib, topN=50):
        sort_pred = np.argsort(prediction)[::-1]
        predicted_templates = {}
        for i in range(1, topN + 1):
            pred_temp = lib[sort_pred[-1][-i]]
            predicted_templates[i] = {'template': pred_temp[-2], 'classification': pred_temp[-3], 'ID': pred_temp[0]}

        return predicted_templates

    def _filter_by_cumulative_predictions(self, actions):
        filtered_actions = []
        for a in actions:
            cumulative = a["cumulative_probability"]
            prob = a['prior']
            if cumulative - prob <= self.cutoff_cumulative:
                filtered_actions.append(a)
        return filtered_actions

    def _get_actions(self, smi):

        cutoff = self.cutoff_number

        prediction = self._get_prediction(smi)
        predicted_templates = self._get_templates(prediction, self.lib, cutoff + 1)
        sort_pred = np.sort(prediction)[::-1]

        actions = []

        for i in range(1, cutoff + 1):
            template = predicted_templates[i]['template']


            metadata = {"policy_probability": round(float(sort_pred[-1][-i]),5),
                        "template_code": predicted_templates[i]['ID'],
                        'classification': predicted_templates[i]['classification']}

            reaction = {'smarts': template,
                        'metadata': metadata,
                        'prior': sort_pred[-1][-i],
                        "cumulative_probability": sum(sort_pred[-1][-i:])}

            actions.append(reaction)

        actions = self._filter_by_cumulative_predictions(actions)

        return actions

    def get_rxns(self, smile):
        if self.policy_model == None:
            self.load_model()

        reactions = self._get_actions(smile)
        rxns = {}
        metadata = {}

        for reaction in reactions:
            name = f"Chem_RB_{reaction['metadata']['classification']}"
            num = 1
            extra_string = f"__{num}"
            while name+extra_string in rxns:
                extra_string = f"__{num}"
                num += 1
            name = name+extra_string
            rxns[name] = [reaction['smarts']]
            metadata[name] = reaction['metadata']
        return rxns, metadata


if __name__ == '__main__':
    action_getter = RingBreaker_ActionGetter()

    target = "c1ccc2c3c([nH]c2c1)CCCCCC3"

    rxns, metadata = action_getter.get_rxns(target)
    print(rxns)
    print(metadata)





