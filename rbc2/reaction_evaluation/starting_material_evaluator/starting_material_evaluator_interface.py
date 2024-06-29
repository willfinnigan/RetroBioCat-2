from abc import ABC, abstractmethod
from typing import Tuple

from rbc2.configs.source_mol_config import SourceMol_Config


class StartingMaterialEvaluatorInterface(ABC):

    def __init__(self, config: SourceMol_Config):
        self.config = config

    @abstractmethod
    def eval(self, smi: str) -> Tuple[bool, dict]:
        """Evaulate the smiles for whether its available as a starting material.
        Returns a tuple of (bool, dict) where the bool is True if the smiles is available and the dict contains
        additional information about the availability.
        """
        pass

    @abstractmethod
    def is_mol_chiral(self, smi: str) -> bool:
        """Return True if the molecule is chiral"""
        pass



