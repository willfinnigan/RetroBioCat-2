from dataclasses import asdict, dataclass
from typing import List, Optional

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.default_expander import DefaultExpander
from rbc2.expansion.expanders.retrobiocat_reaction_retrieval.yaml_rxn_class import YAML_RetroBioCatReactions
from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentData
from rbc2.precedent_identification.data_retrieval.retrobiocat.local_data_query import RetroBioCatLocalPrecedentData
from rbc2.precedent_identification.data_retrieval.retrobiocat.rank_precedents import get_best_enzymes
from rbc2.precedent_identification.similarity_scorer import SimilarityScorer
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluatorInterface
from rbc2.reaction_network_entities.network import Network
from rbc2.reaction_network_entities.reaction import Reaction, sort_reactions_by_score
from rbc2.reaction_network_entities.reaction_option import ReactionOption, sort_options_by_score
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions
from rbc2.template_application.create_reactions_from_output.process_reactions import reaction_sanity_check, \
    process_reactions


@dataclass
class RetroBioCatExpanderConfig():
    enable_precedent_search = True
    search_literature_precedent = True

    include_experimental = False
    include_two_step = True
    include_requires_absence_of_water = False

    score_similarity_before_option_creation = True
    only_active_literature_precedent = True
    similarity_cutoff = 0.55


class RetroBioCatExpander(DefaultExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 starting_material_evaluator: Optional[StartingMaterialEvaluatorInterface] = None,
                 rbc_rxn_class: Optional = None,
                 precedent_data: Optional[PrecedentData] = None,
                 include_experimental: bool = False,
                 include_two_step: bool = True,
                 include_requires_absence_of_water: bool = False,
                 score_similarity_before_option_creation: bool = True,
                 search_literature_precedent: bool = True,
                 only_active_literature_precedent: bool = True,
                 similarity_cutoff: bool = 0.55):

        super().__init__(network=network, config=config)

        self.rxn_type = 'retrobiocat'
        self.rxn_domain = 'biocatalysis'
        self.score_key = ''

        self.rxn_class = rbc_rxn_class
        if self.rxn_class is None:
            self.rxn_class = YAML_RetroBioCatReactions(include_experimental=include_experimental,
                                                       include_two_step=include_two_step,
                                                       include_requires_absence_of_water=include_requires_absence_of_water)

        if precedent_data is None:
            precedent_data = RetroBioCatLocalPrecedentData()

        self.precedent_scorer = SimilarityScorer(precedent_data,
                                                 ranking_function=get_best_enzymes)
        self.starting_material_evaluator = starting_material_evaluator   # if not None, will use substrate availability as part of the option score

        self.score_similarity_before_option_creation = score_similarity_before_option_creation
        self.search_literature_precedent = search_literature_precedent
        self.only_active_literature_precedent = only_active_literature_precedent
        self.similarity_cutoff = similarity_cutoff

    def precedent_evaluation_function(self, reaction: Reaction):
        if self.search_literature_precedent is False:
            return None

        reaction.precedents = self.precedent_scorer.get_precedents(target_smi=reaction.product,
                                                                   topn=1,
                                                                   cutoff=self.similarity_cutoff,
                                                                   reaction_name=reaction.name)

    def create_option(self, smi: str, name: str, smarts: List,
                      template_metadata: dict, score: float) -> ReactionOption:

        option = ReactionOption(smi,
                                name,
                                [],  # empty smarts dict
                                template_metadata,
                                self.rxn_type,
                                self.rxn_domain,
                                score,
                                None,  # none evaluation method because already evaluated
                                precedent_search_function=self.precedent_evaluation_function)
        return option


    def _create_new_options(self, smi: str) -> List[ReactionOption]:
        """Overwrites the default expander method"""

        self.calls += 1
        reactions = self._get_reactions(smi)
        if self.network is not None:
            [self.network.add_reaction(reaction) for reaction in reactions]
        options = []

        search_precedents = (self.score_similarity_before_option_creation == True) and (self.search_literature_precedent == True)

        for reaction in reactions:
            if search_precedents == True:
                self.precedent_evaluation_function(reaction)

            self._calc_reaction_score(reaction)

            # set selected enzyme
            if len(reaction.precedents) > 0:
                reaction.template_metadata[reaction.name]['selected_enzyme'] = reaction.precedents[0].data['enzyme_type']

        reactions = sort_reactions_by_score(reactions)

        if self.config.max_reactions is not None:
            reactions = reactions[:self.config.max_reactions]

        for reaction in reactions:

            option = self.create_option(smi=smi,
                                        name=reaction.name,
                                        smarts=[],
                                        template_metadata=reaction.template_metadata,
                                        score=reaction.score)

            option.evaluated = True
            option.precedents_searched = search_precedents
            option.reactions = [reaction]
            options.append(option)

        sorted_options = sort_options_by_score(options)
        if self.network is not None:
            self.network.bulk_add_options(smi, self.rxn_type, sorted_options)

        return sorted_options

    def _calc_reaction_score(self, reaction: Reaction):
        """Score this option using specificity and complexity"""

        # calculate score
        score = reaction.get_complexity_change() + reaction.get_similarity_score()

        # force between [0, 1] - should we use a flattening function here?
        if score > 1:
            score = 1
        if score < 0:
            score = 0

        # if a biocatalytic step results in all substrates being available, we want to prioritise this
        if isinstance(self.starting_material_evaluator, StartingMaterialEvaluatorInterface):
            if self._are_substrates_all_available(reaction) is True:
                score = 1.01

        # set the reaction score
        reaction.score = score

    def _are_substrates_all_available(self, reaction: Reaction):
        all_available = True
        for smi in reaction.substrates:
            eval, info = self.starting_material_evaluator.eval(smi)
            if eval == 0:
                all_available = False
        return all_available

    def _get_reactions(self, target_smi):
        outcomes = self.rule_applicator.apply_rdchiral(target_smi, self.rxn_class.get_rxns(), multistep_rxns=self.rxn_class.get_multistep_rxns())
        metadata = self._get_metadata(list(outcomes.keys()))
        reactions = create_reactions(target_smi,
                                     outcomes,
                                     metadata=metadata,
                                     rxn_type=self.rxn_type,
                                     rxn_domain=self.rxn_domain)

        reactions = reaction_sanity_check(reactions)

        if self.network is not None:
            reactions = process_reactions(reactions, self.network, self.config)

        return reactions

    def _get_metadata(self, rxn_names: List[str]):
        metadata = {}

        for name in rxn_names:
            metadata[name] = {}
            possible_enzymes = self.rxn_class.get_possible_enzymes(name)
            metadata[name]['possible_enzymes'] = possible_enzymes
            metadata[name]['enzyme_choices'] = [[enz, f"{enz} - {self.rxn_class.get_enzyme_full_name(enz)}"] for enz in possible_enzymes]
            metadata[name]['selected_enzyme'] = metadata[name]['possible_enzymes'][0]
            metadata[name]['enzyme_cofactors'] = self.rxn_class.get_rxn_cofactors(name)
            metadata[name]['is_enzyme'] = 1
            if metadata[name]['selected_enzyme'] == 'Chemical':
                metadata[name]['is_enzyme'] = 0
        return metadata



if __name__ == '__main__':
    expander = RetroBioCatExpander()
    reactions = expander.get_reactions('c1ccc([C@@H]2CCCCN2)cc1')

    for reaction in reactions:
        print()
        print(asdict(reaction))
        print()
