from typing import Optional, Sequence, List

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expander_interface import Expander
from rbc2.expansion.expanders.biosynthesis_expanders import RetroRulesExpander, BKMSExpander, \
    EnzymeMapExpander
from rbc2.expansion.expanders.chemistry_expanders import AIZynthfinderExpander, AskcosPolicyExpander, \
    RingBreakerPolicyExpander
from rbc2.expansion.expanders.retrobiocat_expander import RetroBioCatExpander
from rbc2.data_model.network import Network

RETROBIOCAT = 'retrobiocat'
AIZYNTHFINDER = 'aizynthfinder'
RINGBREAKER = 'ringbreaker'
RETRORULES = 'retrorules'
ASKCOS = 'askcos'
BKMS = 'bkms'
ENZYMEMAP = 'enzymemap'

expander_repo = {RETROBIOCAT:  RetroBioCatExpander,
                 AIZYNTHFINDER: AIZynthfinderExpander,
                 RINGBREAKER: RingBreakerPolicyExpander,
                 RETRORULES: RetroRulesExpander,
                 ASKCOS: AskcosPolicyExpander,
                 BKMS: BKMSExpander,
                 ENZYMEMAP: EnzymeMapExpander
                 }

def get_expanders(expander_names: Sequence[str],
                  network: Optional[Network] = None,
                  expander_config: Optional[Expansion_Config] = None) -> List[Expander]:
    """ Returns expanders given a list of names """

    # if not expansion config specified, use default
    if expander_config is None:
        expander_config = Expansion_Config()

    # get expanders
    expanders = []
    for name in expander_names:
        if name not in expander_repo:
            raise ValueError(f'Expander {name} not found in repository')
        exp = expander_repo[name](config=expander_config,
                                  network=network)
        expanders.append(exp)

    return expanders

