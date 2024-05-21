from collections.abc import Callable

from rbc2.reaction_evaluation.filter_models.aizynthfinder_filter_model import AIZynthFilterModel
from rbc2.data_model.reaction import Reaction


"""
A filter is a function that takes in a reaction, scores it for 'feasability' 
(with 1 being the most feasible) and returns that score.
"""
Filter = Callable[[Reaction], float]

RETROBIOCAT_FILTER = 'retrobiocat_filter'
AIZYNTHFINDER_FILTER = 'aizynthfinder_feasability_filter'

def aizynthfinder_feasability(reaction: Reaction) -> float:
    if reaction.feasability_filter_scores.get(AIZYNTHFINDER_FILTER) is None:
        score = AIZynthFilterModel().run_model(reaction.product, reaction.substrates)
        score = round(score, 4)
        reaction.feasability_filter_scores[AIZYNTHFINDER_FILTER] = score

    return reaction.feasability_filter_scores[AIZYNTHFINDER_FILTER]

def retrobiocat_feasability(reaction: Reaction) -> float:
    if reaction.feasability_filter_scores.get(RETROBIOCAT_FILTER) is None:
        if len(reaction.precedents) == 0:
            score = 0
        else:
            score = max([p.similarity for p in reaction.precedents])
        reaction.feasability_filter_scores[RETROBIOCAT_FILTER] = score

    return reaction.feasability_filter_scores[RETROBIOCAT_FILTER]



default_filter_repo = {AIZYNTHFINDER_FILTER: aizynthfinder_feasability,
                       RETROBIOCAT_FILTER: retrobiocat_feasability}




