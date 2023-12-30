from abc import ABC, abstractmethod
from typing import Tuple

from rbc2.configs.source_mol_config import SourceMol_Config


class StartingMaterialEvaluatorInterface(ABC):

    def __init__(self, config: SourceMol_Config):
        self.config = config

    @abstractmethod
    def eval(self, smi: str) -> Tuple[bool, dict]:
        pass

    @abstractmethod
    def is_mol_chiral(self, smi: str) -> bool:
        pass



