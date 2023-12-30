import pytest

from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction import Reaction


def test_creation_of_one_step_pathway():
    reaction = Reaction('CCCC=O', ['CCCCO'])
    pathway = Pathway([reaction])
    assert pathway.target_smi == 'CCCC=O'

def test_creation_of_two_reaction_pathway_errors_if_different_products():
    reaction1 = Reaction('CCCC=O', ['CCCCO'])
    reaction2 = Reaction('CCCCN', ['CCCCO'])
    with pytest.raises(Exception):
        pathway = Pathway([reaction1, reaction2])

def test_creation_of_two_reaction_pathway_errors_if_cyclic():
    reaction1 = Reaction('CCCC=O', ['CCCCO'])
    reaction2 = Reaction('CCCCO', ['CCCC=O'])
    with pytest.raises(Exception):
        pathway = Pathway([reaction1, reaction2])

def test_pathway_length_of_two_step_pathway_is_2():
    reaction1 = Reaction('CCCC=O', ['CCCCO'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway = Pathway([reaction1, reaction2])
    assert pathway.pathway_length == 2

def test_pathway_length_of_branched_pathway_is_2():
    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    reaction3 = Reaction('C=O', ['CO'])
    pathway = Pathway([reaction1, reaction2, reaction3])
    assert pathway.pathway_length == 2

def test_pathway_end_smi_depths():
    reaction1 = Reaction('CCCC=O', ['CCCCO', 'C=O'])
    reaction2 = Reaction('CCCCO', ['CCCCN'])
    pathway = Pathway([reaction1, reaction2])
    assert pathway.end_smi_depths == {'CCCCN': 2, 'C=O': 1}

def test_target_only_pathway():
    pathway = Pathway([], target_smi='CCCC=O')
    assert pathway.end_smi_depths == {'CCCC=O': 0}



