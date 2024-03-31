from rbc2.expansion.expanders.chemistry_expanders import AIZynthfinderExpander
from rbc2.expansion.expanders.retrobiocat_expander import RetroBioCatExpander
from rbc2.data_model.network import Network


def test_save_load_network():
    network = Network()
    expander = RetroBioCatExpander(network=network)
    reactions = expander.get_reactions('CCCC=O')

    data = network.save()

    network2 = Network()
    expander2 = RetroBioCatExpander(network=network2)
    network2.load(data, [expander2])

    print(network2.reaction_options)

    assert len(network2.all_smis()) == len(network.all_smis())
    assert len(network2.all_reactions()) == len(network.all_reactions())
    assert len(network2.all_reaction_options()) == len(network.all_reaction_options())

def test_loaded_network_can_have_options_evaluated():
    network = Network()
    expander = AIZynthfinderExpander(network=network)
    options = expander.get_options('CCCC=O')

    data = network.save()

    network2 = Network()
    expander2 = AIZynthfinderExpander(network=network2)
    network2.load(data, [expander2])
    options = network2.get_reaction_options('CCCC=O', 'aizynthfinder')
    options[0].evaluate()
    options[1].evaluate()
    assert len(network2.all_reactions()) == 2


