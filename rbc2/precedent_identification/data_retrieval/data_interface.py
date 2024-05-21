from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Tuple
import pandas as pd


class PrecedentData(ABC):

    @abstractmethod
    def __init__(self):
        self.product_column: str = ''
        self.substrate_columns: List[str] = []
        self.enzyme_column: str = ''
        self.id_column: str = ''

    def get_rxn_smi(self, row):

        substrates = []
        for col in self.substrate_columns:
            if col not in list(row.index):  # if column not in row (which is a Series), skip
                continue

            if row[col] is not None and row[col] != '':  # if column is not empty, add to substrates
                substrates.append(row[col])
        substrates = [s for s in substrates if pd.isna(s)==False]

        rxn_smi = f"{'.'.join(substrates)}>>{row[self.product_column]}"  # compose rxn_smi string

        return rxn_smi


    @abstractmethod
    def query_data(self, **key_atts) -> pd.DataFrame():
        """ Takes a number of keyword arguments, specific for each data source"""
        pass

    @abstractmethod
    def get_fps(self, smis: List[str]) -> Tuple[List, List[str]]:
        """ Takes a list of smiles, and returns a list of fingerprints and a list of smiles that were converted to fingerprints """
        pass


