from __future__ import annotations
from dataclasses import asdict
from typing import TYPE_CHECKING

from rbc2.utils.add_logger import add_logger


if TYPE_CHECKING:
    from rbc2.expansion.default_expander_interface import Expander
    from rbc2.reaction_network_entities.reaction_option import ReactionOption

logger = add_logger("SavingLoading_ReactionOptions", level='DEBUG')

def option_to_dict(option: ReactionOption) -> dict:
    """
    Saves the data of a ReactionOption to a dict
    Use the option_from_dict function to recreate the ReactionOption
    """
    opt_dict = {'target_smi': option.target_smi,
                'name': option.name,
                'smarts': option.smarts,
                'metadata': option.metadata,
                'rxn_type': option.rxn_type,
                'rxn_domain': option.rxn_domain,
                'score': option.score,
                'evaluated': option.evaluated,
                'precedents_searched': option.precedents_searched,
                'reaction_ids': [reaction.unique_id for reaction in option.reactions]
                }
    return opt_dict

def option_from_dict(option_dict: dict, expander: Expander) -> ReactionOption:
    """
    Load a previously saved option from a dict - using the expander method
    """

    if option_dict.get('rxn_type', '') != expander.rxn_type:
        logger.warning('Expander rxn_type does not match option_dict rxn_type')

    return expander.create_option(smi=option_dict['target_smi'],
                                  name=option_dict['name'],
                                  smarts=option_dict['smarts'],
                                  template_metadata=option_dict['metadata'],
                                  score=option_dict['score'])