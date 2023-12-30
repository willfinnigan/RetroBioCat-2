from typing import Optional

import numpy as np
from rbc2.utils.add_logger import add_logger
from rbc2.configs.logging_config import logging_config
from rbc2.mcts.tree_node import  MCTS_Node

selection_logger = add_logger('Selection', level=logging_config.mcts_selection)

class Selection():

    def __init__(self):
        self.metrics = {}

    def select(self, node: MCTS_Node, exploration: float) -> Optional[MCTS_Node]:
        node = selection(node, exploration)
        if node is not None:
            if node.option is not None:
                if node.option.rxn_type not in self.metrics:
                    self.metrics[node.option.rxn_type] = 0
                self.metrics[node.option.rxn_type] += 1
        return node

def selection(node: MCTS_Node, exploration: float) -> Optional[MCTS_Node]:
    """ Move through search tree selecting nodes with highest ucb until finding a node with 0 children (ie a leaf)"""

    start_depth = node.depth
    selection_logger.debug(f'-- Selection from depth {start_depth}')

    if node.is_root == True and (node.fully_searched == True or node.terminal == True):
        selection_logger.warning("Search ended at selection because root is fully searched / terminal")
        return None  # if no expansions possible from the root, return None to end the search

    while _should_continue_selection(node):  # loop until node has no children
        node = find_most_promising_child(node, exploration)  # will go to parent if no children

    return node

def _should_continue_selection(node: MCTS_Node):
    if node is None:
        selection_logger.debug(f'-- Selection ended because it returned None')
        return False
    if node.is_root is True and node.fully_searched is True:
        selection_logger.debug(f'-- Selection ended because root is now fully searched')
        return False  # ie is root and has no children left to explore, the search is complete
    if len(node.children) == 0:
        selection_logger.debug(f'-- Selection ended at depth {node.depth} because node has no children')
        return False
    return True


def ucb_calc(node_score: float,
             node_visits: int,
             exploration: float,
             parent_visits: int) -> float:

    # wi/ni + c*sqrt(ln(N)/ni)
    return (node_score / node_visits) + (exploration * (np.sqrt(np.log(parent_visits) / node_visits)))


def ucb_node(node: MCTS_Node, exploration: float):
    if node.terminal or node.fully_searched:
        return -np.inf
    return ucb_calc(node.value, node.visits, exploration, node.parent.visits)


def find_most_promising_child(node: MCTS_Node, exploration: float) -> Optional[MCTS_Node]:
    if node.is_root == True:
        selection_logger.debug('Finding most promising child of root')

    if len(node.children) == 0:  # if no children, this node is terminal, then go back to parent
        node.terminal = True
        selection_logger.debug(f'No children, node is terminal, returning to parent')
        return node.parent

    children_ucbs = [ucb_node(child, exploration) for child in node.children]

    if logging_config.mcts_selection == 'DEBUG':
        debug_scores = [(round(score, 2), node.get_last_rxn_type()) for score, node in zip(children_ucbs, node.children)]
        selection_logger.debug(f'Depth {node.depth} children ucb scores and option type: {debug_scores}')

    max_score = max(children_ucbs)

    if max_score == -np.inf:        # ucb scores are -np.inf if terminal or fully searched
        node.fully_searched = True  # if all -np.inf, then this node has been fully searched
        selection_logger.debug(f'All ucb scores are -np.inf so the node is fully searched, returning to parent')
        return node.parent          # .. in which case we want to go back to the parent

    i = children_ucbs.index(max_score)
    return node.children[i]