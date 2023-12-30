from rbc2.configs.expansion_config import Expansion_Config
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator
from rbc2.reaction_network_entities.network import Network
from rbc2.reaction_network_entities.reaction import Reaction
from rbc2.reaction_network_entities.reaction_option import create_evaluate_option_method, ReactionOption


def get_test_evaluation_func(network=None):
    if network is None:
        network = Network()
    rule_applicator = RuleApplicator()
    expansion_config = Expansion_Config()
    evaluation_function = create_evaluate_option_method(rule_applicator, expansion_config, network=network)
    return evaluation_function

def get_test_option(target, smarts, eval_func=None):
    if eval_func is None:
        eval_func = get_test_evaluation_func()

    option = ReactionOption(target,
                            'test_option',
                            [smarts],
                            {},
                            'test_type',
                            'test_domain',
                            1,
                            eval_func)
    return option

# tests
def test_test_option_creation_is_ok():
    target_smi = 'CCCC=O'
    smarts_str = '[#6:1]=[O:2]>>[#6:1]-[OH:2]'
    option = get_test_option(target_smi, smarts_str)
    assert isinstance(option, ReactionOption)

def test_option_evaluation_returns_reaction():
    target_smi = 'CCCC=O'
    smarts_str = '[#6:1]=[O:2]>>[#6:1]-[OH:2]'
    option = get_test_option(target_smi, smarts_str)
    reactions = option.evaluate()
    assert len(reactions) >= 1
    assert isinstance(reactions[0], Reaction)