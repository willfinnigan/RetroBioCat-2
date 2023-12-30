from rbc2.reaction_network_entities.network import Network
from rbc2.reaction_network_entities.reaction import Reaction
from tests.test_reaction_network_entities.test_reaction_option import get_test_evaluation_func, get_test_option


def test_adding_a_reaction_to_network():
    network = Network()
    reaction = Reaction('CCCC=O', ['CCCCO'],
                        unique_id='test_id',
                        name='test_reaction',
                        rxn_type='retrobiocat',
                        rxn_domain='biocatalysis')
    network.add_reaction(reaction)

    assert network.all_reactions() == [reaction]
    assert sorted(network.all_smis()) == sorted(['CCCC=O', 'CCCCO'])
    assert network.get_reactions_which_molecule_is_substrate_of('CCCCO') == {reaction}
    assert network.get_reactions_which_molecule_is_produced_by('CCCC=O') == {reaction}

def test_option_evaluation():
    target_smi = 'CCCC=O'
    smarts_str = '[#6:1]=[O:2]>>[#6:1]-[OH:2]'

    network = Network()
    eval_func = get_test_evaluation_func(network=network)
    option = get_test_option(target_smi, smarts_str, eval_func=eval_func)


