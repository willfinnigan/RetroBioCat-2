import itertools
from typing import Optional, List, Callable

from rbc2.expansion.default_expander_interface import Expander
from rbc2.expansion.expander_repository import get_expanders
from rbc2.reaction_evaluation.molecular_descriptors import get_mw
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.reaction_option import ReactionOption, sort_options_by_score

combination_method = Callable[[List[List[ReactionOption]]], List[ReactionOption]]

def interleave(list_of_options: List[List[ReactionOption]]) -> List[ReactionOption]:
    interleaved = [x for x in itertools.chain(*itertools.zip_longest(*list_of_options)) if x is not None]
    return interleaved

def order_by_score(list_of_options: List[List[ReactionOption]]) -> List[ReactionOption]:
    flattened = [pathway for sublist in list_of_options for pathway in sublist]
    return sorted(flattened, key=lambda x: x.score, reverse=True)

combination_methods: dict[str: combination_method] = {'interleave': interleave,
                                                      'order_by_score': order_by_score}

class MultiExpander:

    def __init__(self,
                 expanders: dict[str: Expander],
                 network: Optional[Network] = None):

        if len(expanders) == 0:
            raise ValueError("No expanders provided")

        self.expanders = expanders

        # check that all expanders have the same config
        expander_configs = [expander.config for expander in expanders.values()]
        if len(set(expander_configs)) != 1:
            raise ValueError("All expanders must have the same config instance")

        # all expanders should have the same config, so just use the first one
        self.expander_config = list(self.expanders.values())[0].config

        # give all expanders the same network
        for expander in self.expanders.values():
            expander.network = network

    def get_options(self, smis_to_expand: List[str], combination_by: str = 'order_by_score') -> List[ReactionOption]:
        """ For multiple smiles, get the options from each expander and combine them using the combination method """
        per_expander_options = []
        for name, expander in self.expanders.items():
            options = []
            for smi in smis_to_expand:
                if self.is_expander_blocked(smi, expander):
                    continue
                options += expander.get_options(smi)

            options = sort_options_by_score(options)
            per_expander_options.append(options)

        combination_method = combination_methods[combination_by]
        options = combination_method(per_expander_options)
        return options

    def get_reactions(self, smis_to_expand: List[str]) -> List[Reaction]:
        options = self.get_options(smis_to_expand)
        reactions = []
        for opt in options:
            new_reactions = opt.evaluate()
            reactions.extend(new_reactions)
        return reactions

    def template_application_counts(self) -> dict:
        """ Return a dictionary of the number of times a template has been applied for each expander """
        counts = {}
        for expander_name, expander in self.expanders.items():
            counts[expander_name] = expander.number_of_rule_applications()
        counts['total'] = sum([x for x in counts.values()])
        return counts

    def expander_calls(self) -> dict:
        """ Return a dictionary of the number of times each expander has been called """
        counts = {}
        for expander_name, expander in self.expanders.items():
            counts[expander_name] = expander.number_of_calls()
        counts['total'] = sum([x for x in counts.values()])
        return counts

    def is_expander_blocked(self, smi: str, expander: Expander) -> bool:
        """ Return a list of blocked expanders """
        if expander.rxn_domain == 'biocatalysis' or expander.rxn_domain == 'biosynthesis':
            if self.expander_config.use_max_mw_for_enzymes is True:
                if get_mw(smi) > self.expander_config.max_mw_to_use_enzymes:
                    return True
        return False


if __name__ == '__main__':
    expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
    multi_expander = MultiExpander(expanders)
    options = multi_expander.get_options(['CC(=O)OCC(COC(=O)C)OC(=O)C'])
    print(options)
