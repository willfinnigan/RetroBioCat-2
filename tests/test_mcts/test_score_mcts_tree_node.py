from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.mcts_loop.score_node import score_pathway, score_node, penalty_fomula
from rbc2.mcts.tree_node import create_node_with_pathway, create_root
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator import \
    DefaultSQLStartingMaterialEvaluator
from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction import Reaction


def test_calc_complexity_penalty_stays_within_bounds():
    max_penalty = -0.1  # the max penalty that will be returned
    no_penalty = 0
    rel_complexity_no_penalty = 0  # the relative complexity at which there is no penalty
    rel_complexity_max_penalty = -1.5  # the relative complexity at which the max penalty is returned

    for rel_com in [-3,-2, -1.5, -1, -0.75, -0.5, -0.25, 0,1,2]:
        penalty = penalty_fomula(rel_com,
                                 max_penalty,
                                 no_penalty,
                                 rel_complexity_no_penalty,
                                 rel_complexity_max_penalty)

        assert penalty <= no_penalty  # weird inequalities because the penalty is negative
        assert penalty >= max_penalty

def test_score_pathway():
    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway = Pathway([reaction1, reaction2])
    mcts_config = MCTS_Config()
    mcts_config.score_mode = 'basic'
    mcts_config.use_pathway_length_score = False
    score = score_pathway(pathway, mcts_config, DefaultSQLStartingMaterialEvaluator())
    assert score == 1

def test_score_node():
    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway = Pathway([reaction1, reaction2])
    root = create_root('CCCC=O')
    node = create_node_with_pathway(root, pathway, 1)
    mcts_config = MCTS_Config()
    mcts_config.score_mode = 'basic'
    mcts_config.use_pathway_length_score = False
    score = score_node(node, mcts_config, DefaultSQLStartingMaterialEvaluator())
    assert score == 1