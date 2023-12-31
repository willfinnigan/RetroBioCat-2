from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tensorflow import keras
from scipy.special import softmax

from rbc2.configs.download_data_files.download_enzymemap import does_enzymemap_exist, download_enzymemap
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.configs.data_path import path_to_data_folder
from rbc2.utils import load_keras_models
try:
    from tensorflow.python.ops.numpy_ops import np_config
    np_config.enable_numpy_behavior()
except:
    pass


import tensorflow as tf
from functools import partial

data_folder = f'{path_to_data_folder}/enzymemap/brenda'

def sparse_categorical_crossentropy_from_logits(labels, logits):
    return tf.keras.losses.sparse_categorical_crossentropy(labels, logits, from_logits=True)

def top_k(k=1):
    partial_fn = partial(tf.keras.metrics.sparse_top_k_categorical_accuracy, k=k)
    partial_fn.__name__ = 'top_{}'.format(k)
    return partial_fn


def build_model(
    input_shape, output_shape, num_hidden, hidden_size,
    activation='relu', output_activation=None, dropout=0.0, clipnorm=None,
    optimizer=None, learning_rate=0.001,
    compile_model=True, loss=None, metrics=None
):
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.Input(input_shape))
    for _ in range(num_hidden):
        model.add(tf.keras.layers.Dense(hidden_size, activation=activation))
        if dropout:
            model.add(tf.keras.layers.Dropout(dropout))
    model.add(tf.keras.layers.Dense(output_shape, activation=output_activation))
    if optimizer is None or optimizer == 'adam':
        optimizer = tf.keras.optimizers.Adam(learning_rate)
    if clipnorm is not None:
        optimizer.clipnorm = clipnorm
    if compile_model:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics
        )
    return model

def relevance(**kwargs):
    loss = sparse_categorical_crossentropy_from_logits
    metrics = [
        top_k(k=1),
        top_k(k=3),
        top_k(k=5),
        top_k(k=10),
        top_k(k=50),
    ]
    options = {
        'loss': loss,
        'metrics': metrics
    }
    options.update(kwargs)
    return build_model(**options)


class EnzymeMap_Action_Getter():

    def __init__(self,
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 allow_multi_product_templates=False,
                 log_level='WARNING'):
        self.logger = add_logger('EnzymeMapActionGetter', level=log_level)

        if does_enzymemap_exist() == False:
            download_enzymemap()

        self.fp_length = 2048
        self.fp_radius = 2
        self.policy_model = None
        self.templates = None
        self.templates_loaded = False

        self.cutoff_cumulative = cutoff_cumulative
        self.cutoff_number = cutoff_number
        self.allow_multi_product_templates = allow_multi_product_templates

    def load_model(self):
        if self.templates_loaded ==True:
            return

        # Hyperparameters:
        fp_length = 2048
        num_hidden = 1
        hidden_size = 2048
        dropout = 0.2
        learning_rate = 0.001
        activation = 'relu'
        clipnorm = None
        num_classes = 2838

        self.policy_model = relevance(input_shape=(fp_length,),
                          output_shape=num_classes,
                          num_hidden=num_hidden,
                          hidden_size=hidden_size,
                          dropout=dropout,
                          learning_rate=learning_rate,
                          activation=activation,
                          clipnorm=clipnorm
                          )


        policy_path = data_folder + '/default_retro-weights.hdf5'
        self.policy_model.load_weights(policy_path)
        #cls.policy_model = keras.models.load_model(policy_path)  # load_keras_models.LocalKerasModel(policy_path)

        templates_path = data_folder + '/brenda_corrected_template_default_unique_templates.txt'
        with open(templates_path, "r") as file:
            templates_list = file.readlines()
            self.templates = pd.DataFrame(templates_list, columns=['reaction_smarts'])

        self.templates_loaded = True

    def get_rxns(self, smiles):
        self.load_model()

        fp = self._smiles_to_fp(smiles).reshape(1, -1)
        indices, scores = self._predictions_from_model(fp)
        indices, scores = self._filter_by_cumulative_proba(indices, scores)
        indices, scores = self._filter_by_max_number(indices, scores)
        possible_moves = self.templates.iloc[indices]

        smarts_dict, metadata_dict = {}, {}
        for i, (index, row) in enumerate(possible_moves.iterrows()):
            name = f"EM_{i+1}"
            smarts = row['reaction_smarts']
            if self._does_template_have_multiple_products(smarts):
                continue
            metadata = dict(row)
            metadata['score'] = float(scores[i])
            metadata['template_id'] = str(indices[i])
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
    getter = EnzymeMap_Action_Getter()
    smarts_dict, metadata_dict = getter.get_rxns('CCCCC=O')
    print(smarts_dict)
    print(metadata_dict)