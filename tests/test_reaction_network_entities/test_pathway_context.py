from rbc2.pathway_tools.pathway_context_checking import is_previous_reaction_enzyme, \
    count_continuous_previous_chemistry_reactions
from rbc2.reaction_network_entities.pathway import Pathway
from rbc2.reaction_network_entities.reaction import Reaction


def test_previous_reaction_is_enzyme():
    reaction1 = Reaction('CCCCO', ['CCCC=O'], rxn_domain='biocatalysis')
    pathway = Pathway([reaction1])
    assert is_previous_reaction_enzyme(pathway, 'CCCC=O') == True

def test_previous_reaction_is_not_enzyme():
    reaction1 = Reaction('CCCCO', ['CCCC=O'], rxn_domain='chemistry')
    pathway = Pathway([reaction1])
    assert is_previous_reaction_enzyme(pathway, 'CCCC=O') == False


def test_count_continuous_previous_chemistry_reactions_returns_0_for_no_reactions():
    pathway = Pathway([], target_smi='CCCCO')
    assert count_continuous_previous_chemistry_reactions(pathway, 'CCCCO') == 0

def test_count_continuous_previous_chemistry_reactions_returns_1_for_1_reactions():
    reaction1 = Reaction('CCCC=O', ['CCCCO'], rxn_domain='chemistry')
    pathway = Pathway([reaction1])
    assert count_continuous_previous_chemistry_reactions(pathway, 'CCCCO') == 1

def test_count_continuous_previous_chemistry_reactions_returns_2_for_2_reactions():
    reaction1 = Reaction('CCCC=O', ['CCCCO'], rxn_domain='chemistry')
    reaction2 = Reaction('CCCCN', ['CCCC=O'], rxn_domain='chemistry')

    # CCCCO >> CCCC=O >> CCCCN
    pathway = Pathway([reaction1, reaction2])
    assert count_continuous_previous_chemistry_reactions(pathway, 'CCCCO') == 2


def test_count_continuous_previous_chemistry_reactions_returns_1_when_interupted():
    reaction1 = Reaction('CCCC=O', ['CCCCO'], rxn_domain='chemistry')
    reaction2 = Reaction('CCCCN', ['CCCC=O'], rxn_domain='biocatalysis')
    reaction3 = Reaction('CCCC', ['CCCCN'], rxn_domain='chemistry')

    # CCCCO >> CCCC=O >> CCCCN >> CCCC
    pathway = Pathway([reaction1, reaction2, reaction3])
    assert count_continuous_previous_chemistry_reactions(pathway, 'CCCCO') == 1