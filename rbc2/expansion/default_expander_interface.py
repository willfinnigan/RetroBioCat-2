from rbc2.reaction_network_entities.reaction_option import create_evaluate_option_method, sort_options_by_score
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator

from abc import ABC, abstractmethod
from typing import Optional, List
from rbc2.reaction_network_entities.network import Network
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.reaction_network_entities.reaction_option import ReactionOption
from rbc2.reaction_network_entities.reaction import Reaction


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


class DefaultExpander(Expander):
    """A default base class for expanders.  Contains methods which should be overwritten by subclasses,
    and others that optionally can be overwritten"""

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None):

        self.network = network

        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.rule_applicator = RuleApplicator(config=self.config)
        self.option_evaluation_method = create_evaluate_option_method(self.rule_applicator,
                                                                      self.config,
                                                                      self.network,
                                                                      reaction_processing_function=self.reaction_processing_function)

        self.calls = 0

        # these should change
        self.action_getter = None
        self.rxn_type = ''
        self.rxn_domain = ''
        self.score_key = ''

    def precedent_evaluation_function(self, reaction: Reaction) -> bool:
        """ Overwrite this function to evaluate precedents for a reaction"""
        return True

    def reaction_processing_function(self, reactions: List[Reaction]) -> List[Reaction]:
        """ Overwrite this function to apply processing once reactions have been constructed.  Eg to remove cofactors"""
        return reactions

    def get_options(self, smi: str) -> List[ReactionOption]:
        if self.network is None:
            return self._create_new_options(smi)

        if self.network.are_options_available(smi, self.rxn_type):
            options = self.network.get_reaction_options(smi, self.rxn_type)
            sorted_options = sort_options_by_score(options)
            return sorted_options
        else:
            return self._create_new_options(smi)

    def create_option(self,
                      smi,
                      name,
                      smarts,
                      template_metadata,
                      score) -> ReactionOption:

        option = ReactionOption(smi,
                                name,
                                smarts,
                                template_metadata,
                                self.rxn_type,
                                self.rxn_domain,
                                score,
                                self.option_evaluation_method,
                                precedent_search_function=self.precedent_evaluation_function)
        return option

    def get_reactions(self, smi: str) -> List[Reaction]:

        # expand is generate options, then create reactions.
        options = self.get_options(smi)

        reactions = []
        for opt in options:
            new_reactions = opt.evaluate()
            reactions.extend(new_reactions)

        if self.config.max_reactions is not None:
            reactions = reactions[:self.config.max_reactions]

        return reactions

    def number_of_rule_applications(self) -> int:
        return self.rule_applicator.rule_applications

    def number_of_calls(self) -> int:
        return self.calls

    def _create_new_options(self, smi: str):
        self.calls += 1
        smarts_dict, metadata_dict = self.action_getter.get_rxns(smi)
        options = []
        for i, name in enumerate(smarts_dict):
            score = metadata_dict[name][self.score_key]

            option = self.create_option(smi=smi,
                                        name=name,
                                        smarts=smarts_dict[name],
                                        template_metadata={name: metadata_dict[name]},
                                        score=score)

            options.append(option)

        sorted_options = sort_options_by_score(options)

        if self.network is not None:
            self.network.bulk_add_options(smi, self.rxn_type, sorted_options)

        return sorted_options

