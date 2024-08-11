

"""Using the retrobiocat rules sometimes results in pathways with a
unecessary chemical reaction as the first step.  We should filter these out.
"""
from typing import List

from rbc2 import Pathway
from rbc2.pathway_tools.pathway_context_checking import is_enzyme_reaction
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluator


def has_disallowed_first_step_retrobiocat_chemical_reaction(pathway: Pathway, sme: StartingMaterialEvaluator) -> bool:
    """Check if the first step of a pathway is a retrobiocat chemical reaction"""

    # first reaction(s)
    reactions = []
    for smi in pathway.end_smi_depths.keys():
        reactions += list(pathway.smi_substrate_of[smi])

    # are any of these reactions retrobiocat not enzymatic?
    for reaction in reactions:

        # allow all enzyme reactions
        if is_enzyme_reaction(reaction) == True:
            continue
        else:
            # are all the substrates not buyable
            if any([sme.eval(smi)[0] == False for smi in reaction.substrates]):
                return True
            # if the substrates are buyable, is the product also buyable?
            if sme.eval(reaction.product)[0] == True:
                return True
            else:
                continue

    return False

def filter_pathways(pathways: List[Pathway],
                    sme: StartingMaterialEvaluator,
                    verbose=False):

    num_before = len(pathways)
    pathways = [pathway for pathway in pathways if has_disallowed_first_step_retrobiocat_chemical_reaction(pathway, sme) == False]
    num_after = len(pathways)

    if verbose == True:
        print(f"Filtered {num_before-num_after} pathways")
    return pathways