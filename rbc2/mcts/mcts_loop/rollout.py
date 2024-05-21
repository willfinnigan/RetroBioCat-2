from typing import Optional

from rbc2.reaction_evaluation.feasability import Filter
from rbc2.utils.add_logger import add_logger
from rbc2.configs.logging_config import logging_config
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.mcts_loop.evaluate_node import resolve_unevaluated_mcts_node
from rbc2.mcts.mcts_loop.expansion.expand import Expansion
from rbc2.mcts.mcts_loop.selection import Selection
from rbc2.mcts.tree_node import MCTS_Node
from rbc2.data_model.network import Network

rollout_logger = add_logger('Rollout', level=logging_config.mcts_rollout)


def rollout(node: MCTS_Node,
            expansion: Expansion,
            selection: Selection,
            network: Network,
            filters: dict[str: Filter],
            mcts_config: MCTS_Config) -> Optional[MCTS_Node]:
    """
        1. is node terminal
        2. if not, expand the node
        3. selection to go to next node
        4. if node needs evaluating then do this
        5. repeat
    """

    if node is None:
        rollout_logger.debug(f'No rollout because node is None')
        return None

    start_depth = node.depth

    while node.terminal is False and node.fully_searched is False:
        if node.is_evaluated() is False:
            rollout_logger.debug(f'Evaluating node at depth {node.depth}')
            node = resolve_unevaluated_mcts_node(node, network, filters, mcts_config)

        if node.expanded is False:
            expansion.expand(node)

        node = selection.select(node, mcts_config.exploration)

        if node is None:
            return None

    rollout_logger.debug(f'Rollout from depth {start_depth} to depth {node.depth}')
    return node







