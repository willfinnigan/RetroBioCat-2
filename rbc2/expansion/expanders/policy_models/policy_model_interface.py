from abc import ABC, abstractmethod
from typing import List, Tuple

smarts_dict = dict[str: List[str]]
metadata_dict = dict[str: dict]


class PolicyModel(ABC):

    @abstractmethod
    def get_rxns(self, smi: str) -> Tuple[smarts_dict, metadata_dict]:
        pass
