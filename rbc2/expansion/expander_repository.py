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
EGRAPH2RSMILES = 'egraph2rsmiles'
GRAPH2SMILES = 'graph2smiles'
AUGMENTED_TRANSFORMER = 'augmented_transformer'

expander_repo = {RETROBIOCAT:  RetroBioCatExpander,
                 AIZYNTHFINDER: AIZynthfinderExpander,
                 RINGBREAKER: RingBreakerPolicyExpander,
                 RETRORULES: RetroRulesExpander,
                 ASKCOS: AskcosPolicyExpander,
                 BKMS: BKMSExpander,
                 ENZYMEMAP: EnzymeMapExpander,
                 }


def _register_template_free_expanders():
    """Lazily register template-free expanders (require torch)."""
    try:
        from rbc2.expansion.expanders.template_free.template_free_expanders import (
            EGraph2RSmilesExpander, Graph2SmilesExpander, AugmentedTransformerExpander
        )
        expander_repo[EGRAPH2RSMILES] = EGraph2RSmilesExpander
        expander_repo[GRAPH2SMILES] = Graph2SmilesExpander
        expander_repo[AUGMENTED_TRANSFORMER] = AugmentedTransformerExpander
    except ImportError as e:
        import logging
        logging.getLogger(__name__).debug(
            "Template-free expanders not available: %s", e)


_register_template_free_expanders()


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
