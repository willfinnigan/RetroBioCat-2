from typing import List

from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction

def remove_duplicates(reactions: List[Reaction],
                      network: Network,
                      require_same_expander: bool = False,
                      require_same_domain: bool = False,
                      require_matching_name: bool = True,
                      update_match_metadata: bool = False,
                      also_return_matches: bool = False) -> List[Reaction]:

    new_reactions = []

    for reaction in reactions:
        # get reactions with the same product
        same_product_reactions = network.get_reactions_which_molecule_is_produced_by(reaction.product)

        # get reactions which also have the same substrates
        matching_reactions = reactions_which_have_same_substrates(reaction, list(same_product_reactions))

        if require_same_expander == True:
            matching_reactions = reactions_which_have_same_expander(reaction.rxn_type, matching_reactions)
        if require_same_domain == True:
            matching_reactions = reactions_which_have_same_domain(reaction.rxn_domain, matching_reactions)
        if require_matching_name == True:
            matching_reactions = reactions_which_have_same_name(reaction, matching_reactions)

        if len(matching_reactions) == 0:
            new_reactions.append(reaction)
        elif len(matching_reactions) > 1:
            raise Exception('Multiple matching reactions found when looking for duplicates')
        else:
            reaction_match = matching_reactions[0]
            if update_match_metadata is True:
                _merge_template_metadata(reaction_match, reaction)
            if also_return_matches is True:
                new_reactions.append(reaction_match)

    return new_reactions

def _merge_template_metadata(reaction: Reaction, merge_reaction: Reaction):
    for name, metadata in merge_reaction.template_metadata.items():
        while name in reaction.template_metadata:
            name += '_1'
        reaction.template_metadata[name] = metadata

def reactions_which_have_same_substrates(reaction_query: Reaction, reactions_to_check: List[Reaction]) -> List[Reaction]:
    query_substrate_smis = sorted(reaction_query.substrates)
    matching_reactions = []
    for reaction in reactions_to_check:
        substrates_smis = sorted(reaction.substrates)
        if query_substrate_smis == substrates_smis:
            matching_reactions.append(reaction)
    return matching_reactions


def reactions_which_have_same_name(reaction_query: Reaction, reactions_to_check: List[Reaction]) -> List[Reaction]:
    query_name = reaction_query.name
    matching_reactions = []
    for reaction in reactions_to_check:
        check_name = reaction.name
        if check_name == query_name:
            matching_reactions.append(reaction)
    return matching_reactions


def reactions_which_have_same_expander(rxn_type: str, reactions_to_check: List[Reaction]) -> List[Reaction]:
    matching_reactions = []
    for reaction in reactions_to_check:
        if reaction.rxn_type == rxn_type:
            matching_reactions.append(reaction)
    return matching_reactions

def reactions_which_have_same_domain(domain: str, reactions_to_check: List[Reaction]) -> List[Reaction]:
    matching_reactions = []
    for reaction in reactions_to_check:
        if reaction.rxn_domain == domain:
            matching_reactions.append(reaction)
    return matching_reactions