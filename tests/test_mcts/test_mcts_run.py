import time

from rbc2.configs.download_data_files.download_retrorules import does_retrorules_db_exist
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.expansion.expander_repository import get_expanders
from rbc2.expansion.multi_expander import MultiExpander
from rbc2.mcts.mcts_loop.expansion.expand import Expansion
from rbc2.mcts.mcts_loop.rollout import rollout
from rbc2.mcts.mcts_loop.selection import selection, Selection
from rbc2.mcts.mcts import MCTS
from rbc2.mcts.tree_node import create_root
from rbc2.reaction_evaluation.starting_material_evaluator.commercial_starting_material_evaluator import \
    CommercialSME
from rbc2.data_model.network import Network


def test_first_rollout():
    network = Network()
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders, network=network)
    expansion = Expansion(multi_expander, CommercialSME(), MCTS_Config())
    selection = Selection()
    root = create_root('CCCC=O')
    node = selection.select(root, 2)
    filters = {}
    new_node = rollout(node, expansion, selection, network, {hash(root): root}, filters, MCTS_Config())
    assert new_node.is_evaluated() is True

def test_mcts_single_pass():
    target_smi = 'CCCC=O'
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    mcts = MCTS(target_smi, expanders)
    mcts.do_a_loop()

def test_mcts_single_pass_then_selection():
    target_smi = 'CCCC=O'
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    mcts = MCTS(target_smi, expanders)
    mcts.do_a_loop()
    node = selection(mcts.root, 2)
    assert node.depth == 1

def test_run_mcts():
    target_smi = 'CCCC=O'
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    mcts = MCTS(target_smi, expanders)
    mcts.config.max_search_time = 5
    mcts.run()
    assert len(mcts.get_solved_nodes()) > 0


def test_can_get_routes_beyond_buyables():
    target_smi = "Oc1ccc(C=C)cc1F"
    expanders = get_expanders(('retrobiocat',))
    mcts = MCTS(target_smi, expanders)
    mcts.config.max_length = 2
    mcts.config.max_search_time = 5
    mcts.config.allow_moves_beyond_solved = 2
    mcts.run()

    pathways = mcts.get_solved_pathways()
    end_smis = [sorted(pathway.end_smis()) for pathway in pathways]
    assert ['CC(=O)C(=O)O', 'Oc1ccccc1F'] in end_smis


def test_weird_retrorules_case_with_options_but_not_reactions_from_route():
    # mcts should stop once the options are evaluated to nothing

    # only run if we have the necessary retorules files already..
    if does_retrorules_db_exist() == False:
        return

    target_smi = 'CC(O)[C@@H](O)c1ccccc1'
    expansion_config = Expansion_Config()
    expansion_config.rr_diameter = 8
    expansion_config.rr_threshold = 0.4
    expansion_config.rr_max_reactions = 100
    expanders = get_expanders(('retrorules',))
    mcts = MCTS(target_smi, expanders)
    mcts.config.max_length = 4
    mcts.config.max_search_time = 20

    t0 = time.time()
    mcts.run()  # this should be quick because retrorules options evaluate to nothing, so root becomes terminal
    t1 = time.time()
    run_time = t1-t0
    assert run_time < 10












