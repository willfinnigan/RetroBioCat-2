import time
from dataclasses import asdict

import pytest

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expanders.biosynthesis_expanders import BKMSExpander, RetroRulesExpander, EnzymeMapExpander
from rbc2.expansion.expanders.retrobiocat_expander import RetroBioCatExpander
from rbc2.expansion.expanders.chemistry_expanders import AIZynthfinderExpander, AskcosPolicyExpander, RingBreakerPolicyExpander
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator import \
    DefaultSQLStartingMaterialEvaluator
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.reaction_option import ReactionOption

#expanders_to_test = [AIZynthfinderExpander, RingBreakerPolicyExpander, AskcosPolicyExpander, BKMSExpander, RetroRulesExpander, RetroBioCatExpander, EnzymeMapExpander]
expanders_to_test = [RetroBioCatExpander, EnzymeMapExpander, AIZynthfinderExpander]
config = Expansion_Config()
config.rr_diameter = 8
config.rr_combined_score_threshold = 0.3
config.rr_threshold = 0.3
testing_smi = 'C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N'  # tryptophan

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expander_can_create_options(Expander):
    expander = Expander(config=config)
    options = expander.get_options(testing_smi)
    assert isinstance(options[0], ReactionOption)

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expander_options_are_ranked_by_score(Expander):
    expander = Expander(config=config)
    options = expander.get_options(testing_smi)
    for i in range(len(options)-1):
        assert options[i].score >= options[i+1].score

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expander_can_generate_reactions(Expander):
    expander = Expander(config=config)
    reactions = expander.get_reactions(testing_smi)
    assert isinstance(reactions[0], Reaction)

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expander_saves_options_to_network(Expander):
    network = Network()
    expander = Expander(network=network, config=config)
    reactions = expander.get_reactions(testing_smi)
    assert network.are_options_available(testing_smi, expander.rxn_type)

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_options_saved_to_network_are_evaluated(Expander):
    network = Network()
    expander = Expander(network=network, config=config)
    reactions = expander.get_reactions(testing_smi)
    options = network.get_reaction_options(testing_smi, expander.rxn_type)
    assert options[0].evaluated is True

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expander_saves_reactions_to_network(Expander):
    network = Network()
    expander = Expander(network=network, config=config)
    reactions = expander.get_reactions(testing_smi)
    assert len(network.all_reactions()) > 0
    assert len(network.all_smis()) > 0

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_repeat_expander_calls_generates_no_more_reactions(Expander):
    network = Network()
    expander = Expander(network=network, config=config)
    expander.get_reactions(testing_smi)
    num_reactions = len(network.all_reactions())

    # one more call should add no additional reactions because they are already in the network
    expander.get_reactions(testing_smi)
    assert len(network.all_reactions()) == num_reactions


@pytest.mark.parametrize('Expander', expanders_to_test)
def test_expansion_is_fast_if_already_expanded(Expander):
    network = Network()
    expander = Expander(network=network, config=config)
    reactions = expander.get_reactions(testing_smi)
    t0 = time.time()
    reactions = expander.get_reactions(testing_smi)
    t1 = time.time()
    assert t1-t0 < 0.05

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_options_which_are_evaluated_as_having_no_reactions_are_removed_from_network(Expander):
    network = Network()
    expander = Expander(network=network, config=config)

    expander.get_options(testing_smi)

    retrieved_options = network.get_reaction_options(testing_smi, expander.rxn_type)

    options_with_reactions = []
    for option in retrieved_options:
        reactions = option.evaluate()
        if len(reactions) != 0:
            options_with_reactions.append(option)

    new_retrieved_options = network.get_reaction_options(testing_smi, expander.rxn_type)
    assert len(options_with_reactions) == len(new_retrieved_options)

@pytest.mark.parametrize('Expander', expanders_to_test)
def test_can_get_dicts_of_reactions(Expander):
    expander = Expander(config=config)
    reactions = expander.get_reactions(testing_smi)
    reaction_dicts = [asdict(r) for r in reactions]
    assert len(reaction_dicts) > 0

def test_retrobiocat_expander_can_score_based_on_molecule_availability():
    starting_material_evaluator = DefaultSQLStartingMaterialEvaluator()
    starting_material_evaluator.custom_smiles = ['C#CC(O)(CO)CO']
    expander = RetroBioCatExpander(starting_material_evaluator=starting_material_evaluator)
    options = expander.get_options('C#C[C@](O)(C=O)CO')
    assert options[0].name == 'Primary alcohol oxidation'

def test_enzymemap_expander_returns_precedents():
    expander = EnzymeMapExpander()
    expander.config.enzymemap_similarity_cutoff = 0
    reactions = expander.get_reactions('CCCCC=O')

    assert len(reactions[0].precedents) > 0

def test_retrobiocat_multistep_rxns_are_applied():
    expander = RetroBioCatExpander()
    reactions = expander.get_reactions('c1ccc([C@@H]2CCCCN2)cc1')
    print([r.name for r in reactions])
    assert 'Amine deracemization' in [r.name for r in reactions]
