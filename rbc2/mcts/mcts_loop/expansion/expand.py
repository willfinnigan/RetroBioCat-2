from __future__ import annotations

from typing import List

from rbc2.expansion.multi_expander import MultiExpander
from rbc2.mcts.mcts_loop.expansion.expansion_context_filtering import apply_criteria_based_filtering
from rbc2.mcts.mcts_loop.expansion.expansion_score_boost import apply_enzyme_cascade_score_boost
from rbc2.utils.add_logger import add_logger
from rbc2.mcts.mcts_logging_config import logging_config
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.tree_node import MCTS_Node, create_node_from_option, create_root
from rbc2.data_model.reaction_option import ReactionOption
from rbc2.pathway_tools.pathway_evaluation import leaf_molecule_availability


expansion_logger = add_logger('Expansion', level=logging_config.mcts_selection)

class Expansion():

    def __init__(self,
                 multi_expander: MultiExpander,
                 starting_material_evaluator: StartingMaterialEvaluator,
                 mcts_config: MCTS_Config
                 ):
        self.multi_expander = multi_expander
        self.starting_material_evaluator = starting_material_evaluator
        self.mcts_config = mcts_config

    def expand(self, node: MCTS_Node) -> List[MCTS_Node]:
        return expand(node, self.multi_expander, self.starting_material_evaluator, self.mcts_config)


def expand(node: MCTS_Node,
           multi_expander: MultiExpander,
           starting_material_evaluator: StartingMaterialEvaluator,
           mcts_config: MCTS_Config) -> List[MCTS_Node]:

    expansion_logger.debug(f"Expanding node at depth {node.depth} with {len(multi_expander.expanders)} expanders")
    node.expanded = True
    options = get_options(node, multi_expander, starting_material_evaluator, mcts_config)
    node.children = [create_node_from_option(node, option) for option in options]

    # if no children available from expansion, set node to terminal
    if len(node.children) == 0:
        node.terminal = True

    # apply enzyme cascade score boost
    [apply_enzyme_cascade_score_boost(child_node, node.pathway, mcts_config) for child_node in node.children]

    expansion_logger.debug(f"Node has {len(node.children)} children")
    return node.children

def get_options(node: MCTS_Node,
                multi_expander: MultiExpander,
                starting_material_evaluator: StartingMaterialEvaluator,
                mcts_config: MCTS_Config) -> List[ReactionOption]:

    available, not_available = leaf_molecule_availability(node.pathway, starting_material_evaluator)

    if len(not_available) == 0:  # sometimes we might want to continue exploring beyond the buyable stopping criteria
        node.solved = True  # this is the only place solved is updated.
        steps_since_solved = get_mcts_steps_since_solved(node)
        if steps_since_solved < mcts_config.allow_moves_beyond_solved:
            not_available = available

    not_available_and_not_max_length = [smi for smi in not_available if node.pathway.end_smi_depths[smi] <= mcts_config.max_length]

    # if there is a non_available leaf at the max length, maybe there is no point continuing to expand this node..
    if mcts_config.stop_expansion_if_nonbuyable_at_max_length == True:
        if len(not_available_and_not_max_length) < len(not_available):
            expansion_logger.debug(f"Stopping expansion because there is a non-buyable at max length")
            return []

    options = multi_expander.get_options(not_available_and_not_max_length, combination_by=mcts_config.option_combination_method)
    options = apply_criteria_based_filtering(options, node.pathway, mcts_config)
    expansion_logger.debug(f"Options are: {[(opt.rxn_type, opt.score) for opt in options]}")
    return options


def get_mcts_steps_since_solved(node: MCTS_Node) -> int:
    """ When moving beyond a buyable node, we want to know how many steps we have taken since the node was solved """

    if node.solved == False:
        raise Exception('Can not get steps since solved because starting node is not solved')

    steps = -1
    while node.solved is True:
        node = node.parent
        steps += 1
        if node.is_root is True:
            break

    return steps



if __name__ == '__main__':
    from rbc2.expansion.expander_repository import get_expanders

    expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
    multi_expander = MultiExpander(expanders)
    starting_material_evaluator = StartingMaterialEvaluator()

    node = create_root('CCC=O')
    new_nodes = expand(node, multi_expander, starting_material_evaluator, MCTS_Config())




