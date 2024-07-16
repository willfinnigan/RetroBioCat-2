from typing import List

import pandas as pd

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

    if pathway.scores['num_enzyme_steps'] == 0:
        pathway.scores['fraction_enzyme_steps_with_precedent'] = 0
    else:
        pathway.scores['fraction_enzyme_steps_with_precedent'] = pathway.scores['num_enzyme_steps_with_precedent'] / pathway.scores['num_enzyme_steps']

    # a 'simple score' available in the original pathway_explorer code
    if len(pathway.reactions) == 0:
        pathway.scores['pathway_explorer_simple_score'] = 0
    else:
        pathway.scores['pathway_explorer_simple_score'] = ((pathway.scores['change_in_complexity']/pathway.scores['num_steps'])
                                                            + (pathway.scores['change_in_complexity']*pathway.scores['num_enzyme_steps_with_precedent'])) / 2


def pathway_explorer_evaluation(pathways: List[Pathway], sme: StartingMaterialEvaluator, weights=None) -> List[Pathway]:
    """Rank pathways according to the old pathway explorer method (not including the diversity)"""

    for pathway in pathways:  # get the scores if they are not already there
        if 'pathway_explorer_simple_score' not in pathway.scores:
            get_pathway_explorer_scores(pathway, sme)

    if weights is None:  # then use the default weights
        weights = {'normalised_num_enzyme_steps': 1,
                   'normalised_change_in_complexity': 1,
                   'fraction_starting_material_available': 1,
                   'fraction_enzyme_steps_with_precedent': 1}

    df = pd.DataFrame([p.scores for p in pathways])

    df['pathway'] = pathways  # for convenience, have pathway objects in the dataframe

    # calculate normalised scores
    correct_neg_change_in_complexity = 0
    while (df['change_in_complexity'].max() + correct_neg_change_in_complexity) <= 0:
        correct_neg_change_in_complexity += 0.05
    df['normalised_change_in_complexity'] = df['change_in_complexity'] / (
                df['change_in_complexity'].max() + correct_neg_change_in_complexity)

    df['rescaled_num_enzyme_steps'] = df['num_enzyme_steps'].max() - df['num_enzyme_steps']
    df['normalised_num_enzyme_steps'] = df['rescaled_num_enzyme_steps'] / df['rescaled_num_enzyme_steps'].max()

    # calculate weighted scores
    scores_to_df = []
    for index, row in df.iterrows():
        score = 0
        for name in weights:
            score += row[name] * weights[name]
        scores_to_df.append(score)
    df['weighted_score'] = scores_to_df

    df.sort_values('weighted_score', ascending=False, inplace=True)

    # add the scores to the pathways
    for index, row in df.iterrows():
        row['pathway'].scores['weighted_score'] = row['weighted_score']
        row['pathway'].scores['normalised_num_enzyme_steps'] = row['normalised_num_enzyme_steps']
        row['pathway'].scores['normalised_change_in_complexity'] = row['normalised_change_in_complexity']

    # return the pathways in the new sorted order
    return list(df['pathway'])

