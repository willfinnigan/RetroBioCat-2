from typing import List

from rbc2 import Pathway
from rbc2.pathway_tools.pathway_context_checking import is_enzyme_reaction
from rbc2.reaction_evaluation.complexity import get_complexity
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluator

def get_pathway_explorer_scores(pathway: Pathway, sme: StartingMaterialEvaluator):
    """Score a pathway based on the change in complexity"""

    end_smis = pathway.end_smis()
    end_complexities = [get_complexity(smi) for smi in end_smis]
    end_complexity = max(end_complexities)  # Max so that the worst case is considered
    start_complexity = get_complexity(pathway.target_smi)
    pathway.scores['change_in_complexity'] = start_complexity - end_complexity

    pathway.scores['num_steps'] = len(pathway.reactions)
    pathway.scores['num_enzyme_steps'] = sum([int(is_enzyme_reaction(reaction)) for reaction in pathway.reactions])
    pathway.scores['num_steps_with_precedent'] = sum([1 for r in pathway.reactions if len(r.precedents) > 0])
    pathway.scores['num_enzyme_steps_with_precedent'] = sum([1 for r in pathway.reactions if len(r.precedents) > 0 and is_enzyme_reaction(r)==True])

    pathway.scores['fraction_starting_material_available'] = sum([sme.eval(smi)[0] for smi in end_smis]) / len(end_smis)

    # a 'simple score' available in the original pathway_explorer code
    pathway.scores['pathway_explorer_simple_score'] = ((pathway.scores['change_in_complexity']/pathway.scores['num_steps'])
                                                        + (pathway.scores['change_in_complexity']*pathway.scores['num_enzyme_steps_with_precedent'])) / 2


def pathway_explorer_evaluation(pathways: List[Pathway]):
    weights = {'Normalised num enzyme steps': 1,
               'Normalised change in complexity': 1,
               'starting_material': 1,
               'postitive_enzyme_steps_score': 1}
    diversity_weight = 1



