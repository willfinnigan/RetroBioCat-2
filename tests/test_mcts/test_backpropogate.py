from rbc2.mcts.mcts_loop.backpropogate import backpropogate
from rbc2.mcts.tree_node import create_node_with_pathway, create_root
from rbc2.data_model.pathway import Pathway
from rbc2.data_model.reaction import Reaction


def test_score_is_backpropogated_to_the_root():
    root = create_root('CCCC=O')

    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    pathway1 = Pathway([reaction1])
    node1 = create_node_with_pathway(root, pathway1, 1)

    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway2 = Pathway([reaction1, reaction2])
    node2 = create_node_with_pathway(node1, pathway2, 1)

    score = 100

    backpropogate(node2, score)

    assert root.value == 100
    assert root.visits == 2
    assert node1.value == 101
    assert node1.visits == 2
    assert node2.value == 101
    assert node2.visits == 2