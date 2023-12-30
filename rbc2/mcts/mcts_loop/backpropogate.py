from typing import List

from rbc2.mcts.tree_node import MCTS_Node


def backpropogate(node: MCTS_Node, score: float) -> List[MCTS_Node]:
    """
    Backpropogate the score up the tree to the root.
    In the process, collect any new solved nodes and return these
    """

    if node is None:
        return []

    new_solved_nodes = []
    while node is not None:
        node.value += score
        node.visits += 1
        if node.visits == 2 and node.solved == True:  # on its first visit (starts with a 1, and we just added 1), if it is solved, add to solved nodes
            new_solved_nodes.append(node)
        node = node.parent

    new_solved_nodes.reverse()  # shorter pathways are solved first
    return new_solved_nodes


