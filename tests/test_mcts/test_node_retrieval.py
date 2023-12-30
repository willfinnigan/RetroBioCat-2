from rbc2.expansion.expander_repository import get_expanders
from rbc2.mcts.mcts import MCTS


def test_all_returned_nodes_are_evaluated():
    target_smi = '[C@H]1(C2=CC=CC=C2)NCCCC1'
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    mcts = MCTS(target_smi, expanders)  #
    mcts.mcts_config.max_length = 3
    mcts.mcts_config.max_search_time = 5
    mcts.run()

    nodes = mcts.get_all_nodes()
    for node in nodes:
        assert node.is_evaluated()
        assert node.pathway is not None

    print(len(nodes))
    print(mcts.get_run_stats())