from typing import List

import numpy as np

from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.tree_node import MCTS_Node
from rbc2.reaction_evaluation.molecular_descriptors import get_mw
from rbc2.reaction_evaluation.complexity import get_complexity
from rbc2.reaction_evaluation.starting_material_evaluator.commercial_starting_material_evaluator import \
    CommercialSME
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluator
from rbc2.data_model.pathway import Pathway
from rbc2.pathway_tools.pathway_evaluation import leaf_molecule_availability


def score_node(node: MCTS_Node,
               mcts_config: MCTS_Config,
               starting_material_evaluator: StartingMaterialEvaluator):
    if node is None:
        return 0
    return score_pathway(node.pathway, mcts_config, starting_material_evaluator)

def score_pathway(pathway: Pathway,
                  mcts_config: MCTS_Config,
                  starting_material_evaluator: StartingMaterialEvaluator):
    """ Score a pathway based on the number of available molecules in the pathway, optionally using complexity if not available """

    available, not_available = leaf_molecule_availability(pathway, starting_material_evaluator)

    if mcts_config.score_mode == 'basic':
        score = len(available) / len(available + not_available)
    elif mcts_config.score_mode == 'complexity_penalty':
        non_buyable_scores = [mcts_config.non_buyable_score] * len(not_available)
        score = len(available) + sum(non_buyable_scores) / len(available + not_available)  # non_buyable scores are negative or zero
    elif mcts_config.score_mode == 'mass_percent':
        score = calculate_mass_pc_of_available_compounds(available, not_available)
    else:
        raise Exception(f'Unknown score_mode: {mcts_config.score_mode}')

    if mcts_config.use_pathway_length_score:
        transforms_score = 1 / (1 + np.exp(-1 * -(pathway.pathway_length - 4)))  # aizythfinder square function for transforms
        score = (0.95 * score) + (0.05 * transforms_score)

    return score

def calc_non_buyable_penalty_scores(target_smi: str,
                                    not_available: List[str],
                                    mcts_config: MCTS_Config) -> List[float]:
    """ Calculate a penalty based on the complexity of the target and the complexity of the not available molecules """

    target_complexity = get_complexity(target_smi)

    penalty_scores = []
    for smi in not_available:
        relative_complexity = target_complexity - get_complexity(smi)   # if rc is positive, this is good - moving towards simpler materials
        penalty = penalty_fomula(relative_complexity,
                                 mcts_config.max_complexity_penalty,
                                 mcts_config.non_buyable_score,
                                 mcts_config.rel_complexity_no_penalty,
                                 mcts_config.rel_complexity_max_penalty)
        penalty_scores.append(penalty)
    return penalty_scores



def penalty_fomula(rc: float,
                   max_penalty: float,
                   no_penalty: float,
                   rc_no_penalty: float,
                   rc_max_penalty: float,
                   steepness=25) -> float:

    if rc >= rc_no_penalty:
        return no_penalty

    if no_penalty == 0:
        no_penalty = 0.0001  # can't be zero

    sigmoid_range = -rc_max_penalty - rc_no_penalty  # the range over which the function operates
    mid_point = (rc_max_penalty - rc_no_penalty) / 2  # mid point of the sigmoid

    slope = mid_point * steepness * (-sigmoid_range ** -2)
    xoffset = mid_point  # the middle of the sigmoid.
    yoffset = max_penalty  # the max penalty
    y_max = no_penalty  # no penalty

    return (y_max - yoffset) / (1 + np.exp(slope * -(rc - xoffset))) + (yoffset)


def calculate_mass_pc_of_available_compounds(available, not_available):
    if len(not_available) == 0 and len(available) != 0:  # just to save on computation
        return 1

    available_masses = [get_mw(smi) for smi in available]
    non_available_masses = [get_mw(smi) for smi in not_available]
    total_mass = sum(available_masses + non_available_masses)
    available_mass_pc = sum(available_masses) / total_mass
    return available_mass_pc

if __name__ == '__main__':

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

        print(f"Relative complexity = {rel_com}, penalty = {penalty}")

    pc = calculate_mass_pc_of_available_compounds(['CCCCCO'], ['CCC=O'])
    print(pc)

    """
    target_complexity = 2

    for c in [1, 1.25, 1.5, 1.75, 2, 2.5,  2.8, 3.45, 3.5, 4]:
        penalty = calc_complexity_penalty(target_complexity, '',
                                          complexity_for_minus_1,
                                          min_penalty_for_non_buyable,
                                          complexity_for_minimum,
                                          smi_complexity=c)
        print(f"{penalty} penalty at complexity {c}")
    """