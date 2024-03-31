from rbc2.template_application.create_reactions_from_output.check_reaction_direction import does_reaction_go_backwards
from rbc2.data_model.network import Network
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions


def make_test_network():
    network = Network()
    rule_applicator_output = {'test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', rule_applicator_output, rxn_type='test', rxn_domain='not real')
    network.add_reaction(reactions[0])
    return network

def test_backwards_reaction_is_found():

    network = make_test_network()
    example_reaction_outcome = {'backwards_test_rxn': [['CC=O']]}
    reactions = create_reactions('CCO', example_reaction_outcome)
    assert does_reaction_go_backwards(reactions[0], network) == True

def test_normal_reaction_is_false_on_backwards_check():

    network = make_test_network()
    example_reaction_outcome = {'forwards_test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome)
    assert does_reaction_go_backwards(reactions[0], network) == False
