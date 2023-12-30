from __future__ import annotations
import uuid
from dataclasses import dataclass, field
from typing import List, Optional
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.reaction_network_entities.reaction import Reaction
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions
from rbc2.template_application.create_reactions_from_output.process_reactions import process_reactions, reaction_sanity_check
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator

from typing import TYPE_CHECKING

from rbc2.utils.add_logger import add_logger

if TYPE_CHECKING:
    from rbc2.reaction_network_entities.network import Network

logger = add_logger('ReactionOptions', level='DEBUG')

@dataclass
class ReactionOption():
    target_smi: str
    name: str
    smarts: list
    metadata: dict
    rxn_type: str
    rxn_domain: str
    score: float

    # functions which are used to evaluate the option
    evaluation_function: callable
    precedent_search_function: Optional[callable] = None

    # reaction outputs from evaluation
    reactions: List[Reaction] = field(default_factory=list)
    evaluated: bool = False
    precedents_searched: bool = False

    def __post_init__(self):
        self.unique_id = f"{self.name}_{str(uuid.uuid4())}"

    def evaluate(self):
        if self.evaluated is False and self.evaluation_function is not None:
            self.reactions = self.evaluation_function(self)
            self.evaluated = True

        if self.precedents_searched is False and self.precedent_search_function is not None:
            for reaction in self.reactions:
                self.precedent_search_function(reaction)
            self.precedents_searched = True

        return self.reactions





def create_evaluate_option_method(rule_applicator: RuleApplicator,
                                  config: Expansion_Config,
                                  network: Optional[Network] = None,
                                  reaction_processing_function: Optional[callable] = None):

    def evaluate_method(option: ReactionOption):
        """ Evaluates the option to return reactions """

        rdchiral_rxns = rule_applicator.smarts_to_rdchiral(option.smarts)
        rxn = {option.name: rdchiral_rxns}
        outcomes = rule_applicator.apply_rdchiral(option.target_smi, rxn)
        reactions = create_reactions(option.target_smi,
                                     outcomes,
                                     score=option.score,
                                     metadata=option.metadata,
                                     rxn_type=option.rxn_type,
                                     rxn_domain=option.rxn_domain)

        if reaction_processing_function is not None:
            reactions = reaction_processing_function(reactions)

        reactions = reaction_sanity_check(reactions)

        if network is not None:
            reactions = process_reactions(reactions, network, config)
            if len(reactions) == 0:
                network.remove_option(option)
            else:
                [network.add_reaction(reaction) for reaction in reactions]

        return reactions

    return evaluate_method


def sort_options_by_score(options: List[ReactionOption]) -> List[ReactionOption]:
    return sorted(options, key=lambda x: x.score, reverse=True)