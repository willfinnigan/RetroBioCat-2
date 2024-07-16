import json
from pathlib import Path

from rbc2 import Network, Pathway, CommercialSME
from rbc2.pathway_tools.pathway_explorer.grouping import group_pathways
from rbc2.pathway_tools.pathway_explorer.scoring import get_pathway_explorer_scores

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

    print(pathway.__hash__())

def test_grouping_pathways():
    pathway = load_test_pathway()
    get_pathway_explorer_scores(pathway, CommercialSME())

    groups = group_pathways([pathway, pathway, pathway])

    assert len(groups) == 1
