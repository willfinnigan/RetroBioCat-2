from rbc2.mcts.tree_node import create_root, MCTS_Node, create_node_from_option, \
    create_node_with_pathway
from rbc2.data_model.pathway import Pathway
from rbc2.data_model.reaction import Reaction
from tests.test_reaction_network_entities.test_reaction_option import get_test_option


def test_can_create_a_root_node():
    target_smi = 'CCCC=O'
    root = create_root(target_smi)
    assert isinstance(root, MCTS_Node)

def test_can_create_a_node_from_a_reaction_option():
    target_smi = 'CCCC=O'
    smarts_str = '[#6:1]=[O:2]>>[#6:1]-[OH:2]'

    option = get_test_option(target_smi, smarts_str)
    root = create_root(target_smi)

    node = create_node_from_option(root, option)

    assert isinstance(node, MCTS_Node)

def test_can_create_a_node_from_pathway():

    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway = Pathway([reaction1, reaction2])

    root = create_root('CCCC=O')

    node = create_node_with_pathway(parent=root, pathway=pathway, value=1)

    assert isinstance(node, MCTS_Node)
