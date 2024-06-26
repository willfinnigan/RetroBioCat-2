from rbc2.configs import Expansion_Config
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.reaction_option import ReactionOption, sort_options_by_score, create_evaluate_option_method
from rbc2.expansion.expander_interface import Expander
from rbc2.expansion.expanders.policy_models.policy_model_interface import PolicyModel
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator

from typing import Optional, List



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
        self.policy_model: Optional[PolicyModel] = None  # update this with a policy model for selecting reactions
        self.rxn_type = ''
        self.rxn_domain = ''
        self.score_key = ''

    def reaction_processing_function(self, reactions: List[Reaction]) -> List[Reaction]:
        """ Overwrite this function to apply processing once reactions have been constructed.  Eg to remove cofactors"""
        return reactions

    def precedent_evaluation_function(self, reaction: Reaction):
        """ Overwrite this function to evaluate precedents for a reaction.
        Should update the reaction.precedents with new precedents"""
        return True

    def final_processing_function(self, reactions: List[Reaction]) -> List[Reaction]:
        """ Overwrite this function to apply processing to the final list of reactions before returning"""
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

        reactions = self.final_processing_function(reactions)

        return reactions

    def number_of_rule_applications(self) -> int:
        return self.rule_applicator.rule_applications

    def number_of_calls(self) -> int:
        return self.calls

    def _create_new_options(self, smi: str):
        self.calls += 1
        smarts_dict, metadata_dict = self.policy_model.get_rxns(smi)
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

