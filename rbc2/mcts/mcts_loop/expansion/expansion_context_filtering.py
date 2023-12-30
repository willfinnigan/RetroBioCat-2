from typing import List, Optional
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.pathway_tools.pathway_context_checking import number_of_chemistry_reactions, \
    is_previous_reaction_chemistry, count_continuous_previous_chemistry_reactions, \
    count_continous_reactions
from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction_option import ReactionOption


def filter_by_max_chemistry_nodes(options: List[ReactionOption],
                                  pathway: Pathway,
                                  max_chemistry: int) -> List[ReactionOption]:
    """ If there is a max chemistry option, ensure we don't go over this"""

    if max_chemistry is None:
        return options

    # if we're at max chemistry nodes, remove all chemistry options
    if number_of_chemistry_reactions(pathway) >= max_chemistry:  # no more chemistry allowed
        options = [opt for opt in options if opt.rxn_domain != 'chemistry']

    return options


def apply_beginning_and_end_chemistry_filter(options: List[ReactionOption],
                                             pathway: Pathway,
                                             max_chemistry_at_end: Optional[int],
                                             max_chemistry_at_beginning: Optional[int]) -> List[ReactionOption]:


    if len(pathway.reactions) == 0:
        return options  # if no reactions yet, allow anything

    allowed_options = []
    for option in options:

        num_rxns = count_continous_reactions(pathway, option.target_smi)
        chem_rxns = count_continuous_previous_chemistry_reactions(pathway, option.target_smi)


        if num_rxns == chem_rxns:
            # if all previous continuous reactions are chemistry,
            # then allow more chemistry only if we are not over the max allowed chemistry at end
            if option.rxn_domain != 'chemistry':
                allowed_options.append(option)
            elif max_chemistry_at_end is None:
                allowed_options.append(option)
            elif chem_rxns < max_chemistry_at_end:
                allowed_options.append(option)


        elif is_previous_reaction_chemistry(pathway, option.target_smi) == False:
            # must be middle if previous reaction was not chemistry, and not all reactions are chemistry
            allowed_options.append(option)

        else:
            # must be beginning if previous reaction was chemistry, and not all reactions are chemistry
            # only chemistry is now allowed
            if option.rxn_domain == 'chemistry':
                if max_chemistry_at_beginning is None:
                    allowed_options.append(option)
                elif chem_rxns < max_chemistry_at_beginning:
                    allowed_options.append(option)

    return allowed_options


def apply_criteria_based_filtering(options: List[ReactionOption],
                                   pathway: Pathway,
                                   mcts_config: MCTS_Config) -> List[ReactionOption]:
    """ Remove any options based on criteria"""

    # if parent_node has no reactions, no filtering necessary
    if len(pathway.reactions) == 0:
        return options

    # overall max chemistry nodes allowed
    options = filter_by_max_chemistry_nodes(options, pathway, mcts_config.max_chemistry_nodes)

    # filtering if there are limitations on use of chemistry and beginning/end
    if mcts_config.chemistry_only_at_beginning_or_end == True:
        options = apply_beginning_and_end_chemistry_filter(options,
                                                           pathway,
                                                           mcts_config.max_chemistry_at_end,
                                                           mcts_config.max_chemistry_at_beginning)
    return options



