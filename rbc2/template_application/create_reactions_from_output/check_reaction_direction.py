from __future__ import annotations
from typing import List

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from rbc2.reaction_network_entities.network import Network
    from rbc2.reaction_network_entities.reaction import Reaction


def does_reaction_go_backwards(reaction: Reaction, network: Network) -> bool:
    """Returns True if reaction goes backwards"""

    smis_to_check = [smi for smi in reaction.substrates if smi in network.all_smis()]
    for smi in smis_to_check:
        reactions_which_produce_smi = network.get_reactions_which_molecule_is_produced_by(smi)
        for producing_reaction in reactions_which_produce_smi:
            if reaction.product in [smi for smi in producing_reaction.substrates]:
                return True
    return False


def remove_backwards_reactions(reactions: List[Reaction], network: Network):
    return [r for r in reactions if does_reaction_go_backwards(r, network) is False]
