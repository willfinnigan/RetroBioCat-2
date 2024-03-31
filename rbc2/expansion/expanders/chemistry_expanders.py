from typing import Optional

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expanders.action_getters.aizynthfinder.aizynthfinder_actions import \
    AizynthfinderActionGetter
from rbc2.expansion.expanders.action_getters.askcos.askcos_action_getter import Askcos_Action_Getter
from rbc2.expansion.expanders.action_getters.ring_breaker.ringbreaker_actions import \
    RingBreaker_ActionGetter
from rbc2.expansion.default_expander_interface import DefaultExpander
from rbc2.data_model.network import Network


class AIZynthfinderExpander(DefaultExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 template_column='retro_template',
                 cutoff_cumulative=0.995,
                 cutoff_number=50):
        super().__init__(network=network, config=config)
        self.action_getter = AizynthfinderActionGetter(template_column=template_column,
                                                       cutoff_cumulative=cutoff_cumulative,
                                                       cutoff_number=cutoff_number)
        self.rxn_type = 'aizynthfinder'
        self.rxn_domain = 'chemistry'
        self.score_key = 'policy_probability'


class RingBreakerPolicyExpander(DefaultExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 cutoff_cumulative=0.995,
                 cutoff_number=10):
        super().__init__(network=network, config=config)

        self.action_getter = RingBreaker_ActionGetter(cutoff_number=cutoff_number,
                                                      cutoff_cumulative=cutoff_cumulative,
                                                      log_level='WARNING')
        self.rxn_type = 'ringbreaker'
        self.rxn_domain = 'chemistry'
        self.score_key = 'policy_probability'


class AskcosPolicyExpander(DefaultExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 allow_multi_product_templates=False):

        super().__init__(network=network, config=config)
        self.action_getter = Askcos_Action_Getter(cutoff_number=cutoff_number,
                                                  cutoff_cumulative=cutoff_cumulative,
                                                  allow_multi_product_templates=allow_multi_product_templates)
        self.rxn_type = 'askcos'
        self.rxn_domain = 'chemistry'
        self.score_key = 'score'


if __name__ == '__main__':
    #expander = AIZynthfinderExpander()
    #reactions = expander.get_reactions('CCCO')
    #print(reactions)

    expander = AskcosPolicyExpander()
    reactions = expander.get_reactions('CCCCO')
    print(reactions)

    # expander = RingBreakerExpander()
    # reactions = expander.get_reactions('OCC1CCCC(=O)O1')
    # print(reactions)
