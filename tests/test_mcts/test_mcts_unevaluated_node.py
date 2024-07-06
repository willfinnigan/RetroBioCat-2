from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.mcts_loop.evaluate_node import resolve_unevaluated_mcts_node
from rbc2.mcts.tree_node import create_root, MCTS_Node
from rbc2.data_model.network import Network
from tests.test_reaction_network_entities.test_reaction_option import get_test_evaluation_func, get_test_option


def test_evaluate_node():
    target_smi = 'CCCC=O'
    smarts_str = '[#6:1]=[O:2]>>[#6:1]-[OH:2]'

    network = Network()
    eval_func = get_test_evaluation_func(network=network)
    option = get_test_option(target_smi, smarts_str, eval_func=eval_func)

    root = create_root(target_smi)
    child_1 = MCTS_Node(parent=root, option=option, value=1)
    root.children = [child_1]

    filters = {}
    new_child = resolve_unevaluated_mcts_node(child_1, network, {}, filters, MCTS_Config())

    assert new_child != child_1
    assert root.children == [new_child]
    assert new_child.is_evaluated() is True


# some more tests of the evaluate methods would be good here


