from abc import ABC, abstractmethod
from typing import List

smarts_list = List[str]
rxns_dict = dict[str: smarts_list]
steps = List[smarts_list]
multi_rxns_dict = dict[str: List[steps]]

class RetroBioCatReactions(ABC):
    """ By default rbc2 will use the yaml_rxn_class.  However, the same interface (specified here),
    can be used to pull reactions from a database - like the mongodb used by retrobiocat.com"""

    def __init(self,
               include_experimental=False,
               include_two_step=True,
               include_requires_absence_of_water=False,
               reverse=False,
               use_rdchiral=True
               ):

        self.include_experimental = include_experimental
        self.include_two_step = include_two_step
        self.include_requires_absence_of_water = include_requires_absence_of_water
        self.reverse = reverse
        self.use_rdchiral = use_rdchiral

    @abstractmethod
    def get_rxns(self) -> rxns_dict:
        pass

    def get_multistep_rxns(self) -> multi_rxns_dict:
        pass

    @abstractmethod
    def get_possible_enzymes(self, rxn_name) -> List[str]:
        pass

    @abstractmethod
    def get_rxn_cofactors(self, rxn_name) -> dict:
        pass

    @abstractmethod
    def get_enzyme_full_name(self, enzyme) -> str:
        pass





