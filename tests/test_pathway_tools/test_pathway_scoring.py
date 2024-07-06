import json
from dataclasses import asdict
from pathlib import Path

from rbc2 import Network, Pathway, CommercialSME
from rbc2.pathway_tools.pathway_context_checking import is_enzyme_reaction
from rbc2.pathway_tools.pathway_explorer_evaluation import get_pathway_explorer_scores

test_data = Path(__file__).parents[1] / 'test_data'

def load_test_pathway():

    # print contents of test_data
    for item in test_data.iterdir():
        print(item)

    filepath = f"{test_data}/piperidine_network.json"
    with open(filepath, 'r') as f:
        network_data = json.load(f)

    network = Network()
    network.load(network_data, [])
    pathway = Pathway(list(network.reactions))
    return pathway

def test_scoring_of_pathway():
    pathway = load_test_pathway()

    get_pathway_explorer_scores(pathway, CommercialSME())

    assert pathway.scores['change_in_complexity'] == 1.098
    assert pathway.scores['num_steps'] == 4
    assert pathway.scores['num_enzyme_steps'] == 3
    assert pathway.scores['num_steps_with_precedent'] == 4
    assert pathway.scores['num_enzyme_steps_with_precedent'] == 3
    assert pathway.scores['fraction_starting_material_available'] == 1.0


