from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Optional, List, TYPE_CHECKING
from rbc2.configs.expansion_config import Expansion_Config

if TYPE_CHECKING:
    from rbc2 import Network, ReactionOption, Reaction

class Expander(ABC):
    """The expander interface, which defines the methods that an expander must implement"""

    @abstractmethod
    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None):
        self.network = network
        self.rxn_type = ''
        self.rxn_domain = ''
        self.config = config

    @abstractmethod
    def get_options(self, smi: str) -> List[ReactionOption]:
        pass

    @abstractmethod
    def create_option(self, smi: str, name: str, smarts: List[str],
                      template_metadata: dict, score: float) -> ReactionOption:
        pass

    @abstractmethod
    def get_reactions(self, smi: str) -> List[Reaction]:
        pass

    @abstractmethod
    def number_of_rule_applications(self) -> int:
        pass

    @abstractmethod
    def number_of_calls(self) -> int:
        pass