from rbc2.configs.mcts_config import MCTS_Config
from rbc2.pathway_tools.pathway_context_checking import is_previous_reaction_enzyme
from rbc2.mcts.tree_node import MCTS_Node
from rbc2.reaction_network_entities.pathway import Pathway


def apply_enzyme_cascade_score_boost(child_node: MCTS_Node, pathway: Pathway, mcts_config: MCTS_Config):
    """ Update the child nodes score if it should be boosted based on the context of the parent node """

    # if no reactions currently, no boost necessary
    if len(pathway.reactions) == 0:
        return

    # if boosting casades is turned off, no boost needed
    if mcts_config.boost_enzyme_score_if_in_cascade == False:
        return

    # if child node is not an enzyme reaction option, no boost needed
    if child_node.option.rxn_domain not in ['biosynthesis', 'biocatalysis']:
        return

    # if the previous reaction was an enzyme reaction, apply boost
    leaf_smi = child_node.option.target_smi
    if is_previous_reaction_enzyme(pathway, leaf_smi) == True:
        child_node.value += mcts_config.boost_enzyme_in_cascade_score_by