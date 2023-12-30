from typing import List

from rbc2.configs.expansion_config import Expansion_Config

from rbc2.reaction_network_entities.reaction import Reaction
from rbc2.template_application.create_reactions_from_output.check_reaction_direction import \
    remove_backwards_reactions
from rbc2.template_application.create_reactions_from_output.find_duplicates import remove_duplicates
from rbc2.reaction_network_entities.network import Network



def reaction_sanity_check(reactions: List[Reaction]) -> List[Reaction]:
    # ensure all reactions have substrates
    reactions = [reaction for reaction in reactions if len(reaction.substrates) != 0]
    return reactions


def process_reactions(reactions: List[Reaction], network: Network, config: Expansion_Config) -> List[Reaction]:
    if config.allow_backwards == False:
        reactions = remove_backwards_reactions(reactions, network)
    if config.allow_duplicates == False:
        reactions = remove_duplicates(reactions,
                                      network,
                                      require_same_expander=config.duplicates_require_same_expander,
                                      require_same_domain=config.duplicates_require_same_domain,
                                      require_matching_name=config.duplicates_require_same_name,
                                      update_match_metadata=config.merge_duplicate_metadata)



    return reactions