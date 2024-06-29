from abc import ABC, abstractmethod
from typing import Tuple

class StartingMaterialEvaluator(ABC):

    def __init__(self):
        self.target_always_not_buyable = True

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



