from rbc2.expansion.multi_expander import MultiExpander

from rbc2.configs.mcts_config import MCTS_Config
from rbc2.expansion.expander_repository import get_expanders
from rbc2.mcts.mcts_loop.expansion.expand import  get_mcts_steps_since_solved, expand
from rbc2.mcts.tree_node import create_root, MCTS_Node, create_node_with_pathway
from rbc2.reaction_evaluation.starting_material_evaluator import StartingMaterialEvaluator
from rbc2.reaction_network_entities.network import Network
from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction import Reaction


def test_mcts_expansion_returns_unevaluated_mcts_nodes():
    root = create_root('CCCC=O')
    network = Network()
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders, network=network)
    starting_material_evaluator = StartingMaterialEvaluator()
    mcts_config = MCTS_Config()

    new_nodes = expand(root, multi_expander, starting_material_evaluator, mcts_config)
    assert isinstance(new_nodes[0], MCTS_Node)


def test_steps_since_solved_returns_num_nodes_since_solved_is_0_if_only_just_solved():
    root = create_root('C(N)C(N)C(N)C=O')
    reaction1 = Reaction('C(N)C(N)C(N)C=O', ['C(N)C(N)C(N)CO'])
    pathway1 = Pathway([reaction1])
    node_1 = create_node_with_pathway(root, pathway1, 1)
    reaction2 = Reaction('C(N)C(N)C(N)CO', ['CCCO'])
    pathway2 = Pathway([reaction1, reaction2])
    node_2 = create_node_with_pathway(node_1, pathway2, 1)
    node_2.solved = True

    steps = get_mcts_steps_since_solved(node_2)

    assert steps == 0

def test_steps_since_solved_returns_num_nodes_since_solved_is_1_if_one_other_solved_node():
    root = create_root('C(N)C(N)C(N)C=O')
    reaction1 = Reaction('C(N)C(N)C(N)C=O', ['C(N)C(N)C(N)CO'])
    pathway1 = Pathway([reaction1])
    node_1 = create_node_with_pathway(root, pathway1, 1)
    reaction2 = Reaction('C(N)C(N)C(N)CO', ['CCCO'])
    pathway2 = Pathway([reaction1, reaction2])
    node_2 = create_node_with_pathway(node_1, pathway2, 1)
    node_2.solved = True
    reaction3 = Reaction('CCCO', ['CCCN'])
    pathway3 = Pathway([reaction1, reaction2, reaction3])
    node_3 = create_node_with_pathway(node_2, pathway3, 1)
    node_3.solved = True

    steps = get_mcts_steps_since_solved(node_3)

    assert steps == 1


def test_0_allowed_moves_beyond_solved_gives_0_new_nodes():
    network = Network()
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders, network=network)
    starting_material_evaluator = StartingMaterialEvaluator()
    mcts_config = MCTS_Config()

    root = create_root('CCCC=O')
    reaction1 = Reaction('CCCC=O', ['CCCCO'])
    pathway = Pathway([reaction1])
    node_1 = create_node_with_pathway(root, pathway, 1)

    new_nodes = expand(node_1, multi_expander, starting_material_evaluator, mcts_config)

    assert len(new_nodes) == 0

def test_allow_moves_beyond_solved_1_gives_new_nodes_even_though_node_is_solved():
    network = Network()
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders, network=network)
    starting_material_evaluator = StartingMaterialEvaluator()
    mcts_config = MCTS_Config()
    mcts_config.allow_moves_beyond_solved = 1

    root = create_root('CCCC=O')
    reaction1 = Reaction('CCCC=O', ['CCCCO'])
    pathway = Pathway([reaction1])
    node_1 = create_node_with_pathway(root, pathway, 1)

    new_nodes = expand(node_1, multi_expander, starting_material_evaluator, mcts_config)
    assert len(new_nodes) > 0
