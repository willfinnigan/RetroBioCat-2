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
        self.rxn_type = ''  # The type of expander, eg 'retrobiocat' or 'aizynthfinder'
        self.rxn_domain = ''  # The domain of this expander, eg 'biocatalysis' or 'chemistry'
        self.config = config

    @abstractmethod
    def get_reactions(self, smi: str) -> List[Reaction]:
        """Return a list of Reaction objects for a given SMILES string.
        This is main method utilised by an expander"""
        pass

    @abstractmethod
    def get_options(self, smi: str) -> List[ReactionOption]:
        """Return a list of ReactionOption objects for a given SMILES string"""
        pass

    @abstractmethod
    def number_of_rule_applications(self) -> int:
        """Get statistics on the use of this expander"""
        pass

    @abstractmethod
    def number_of_calls(self) -> int:
        """Get statistics on the use of this expander"""
        pass

    @abstractmethod
    def create_option(self, smi: str, name: str, smarts: List[str],
                      template_metadata: dict, score: float) -> ReactionOption:
        """Create a ReactionOption object given the necessary information.
        This is primarily used to load a ReactionOption from a database."""
        pass