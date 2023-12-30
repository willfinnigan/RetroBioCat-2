from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction import Reaction


def is_enzyme_reaction(reaction: Reaction) -> bool:
    return reaction.rxn_domain in ['biosynthesis', 'biocatalysis']

def is_chemistry_reaction(reaction: Reaction) -> bool:
    return reaction.rxn_domain == 'chemistry'

def is_previous_reaction_enzyme(pathway: Pathway, leaf_smi: str) -> bool:
    """ Return True if the previous reaction in the pathway is an enzyme reaction"""

    # if pathway has no reactions, return False
    if len(pathway.reactions) == 0:
        return False

    # get preceeding reaction using the leaf smi of interest
    preceeding_reaction = pathway.get_reaction_with_substrate(leaf_smi)

    # if its an enzyme_reaction -> True
    return is_enzyme_reaction(preceeding_reaction)

def is_previous_reaction_chemistry(pathway: Pathway, leaf_smi: str) -> bool:
    """ Return True if the previous reaction in the pathway is a chemistry reaction"""

    # if pathway has no reactions, return False
    if len(pathway.reactions) == 0:
        return False

    # get preceeding reaction using the leaf smi of interest
    preceeding_reaction = pathway.get_reaction_with_substrate(leaf_smi)

    # if its an enzyme_reaction -> True
    return is_chemistry_reaction(preceeding_reaction)

def count_continuous_previous_chemistry_reactions(pathway: Pathway, leaf_smi: str) -> int:
    """ Recursive function to return the number of continuous chemistry reactions in the pathway """

    # get preceeding reaction using the leaf smi of interest
    preceeding_reaction = pathway.get_reaction_with_substrate(leaf_smi)
    if preceeding_reaction is None:
        return 0

    if is_chemistry_reaction(preceeding_reaction):
        return 1 + count_continuous_previous_chemistry_reactions(pathway, preceeding_reaction.product)
    else:
        return 0

def count_continous_reactions(pathway: Pathway, leaf_smi: str) -> int:
    """ Recursive function to return the number of continuous reactions in the pathway from leaf"""

    # get preceeding reaction using the leaf smi of interest
    preceeding_reaction = pathway.get_reaction_with_substrate(leaf_smi)
    if preceeding_reaction is None:
        return 0

    return 1 + count_continous_reactions(pathway, preceeding_reaction.product)

def get_chemistry_reactions(pathway: Pathway) -> list:
    """ Return the chemistry reactions in the pathway """
    return [rxn for rxn in pathway.reactions if rxn.rxn_domain == 'chemistry']

def get_non_chemistry_reactions(pathway: Pathway) -> list:
    """ Return the non-chemistry reactions in the pathway """
    return [rxn for rxn in pathway.reactions if rxn.rxn_domain != 'chemistry']

def number_of_chemistry_reactions(pathway: Pathway) -> int:
    """ Return the number of chemistry reactions in the pathway """
    return len(get_chemistry_reactions(pathway))