from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.mcts_loop.expansion.expansion_context_filtering import filter_by_max_chemistry_nodes, apply_beginning_and_end_chemistry_filter
from rbc2.data_model.pathway import Pathway
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.reaction_option import ReactionOption


def test_three_chemistry_steps_allows_anything():
    pathway = Pathway([Reaction('CCCO', ['CCC=O'], rxn_domain='chemistry'),
                       Reaction('CCC=O', ['CCC(=O)O'], rxn_domain='chemistry'),
                       Reaction('CCC(=O)O', ['CCC(=O)NC'], rxn_domain='chemistry')])

    option1 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'chemistry', 0.0, lambda x: x)
    option2 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'biocatalysis', 0.0, lambda x: x)

    options = [option1, option2]

    # no filtering because only 3 chemistry steps
    assert filter_by_max_chemistry_nodes(options, pathway, 4) == [option1, option2]

    # filtering because max 3
    assert filter_by_max_chemistry_nodes(options, pathway, 3) == [option2]

    # we've only got chemistry, but beginning is over max so can only add biocatalysis
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    2, 2) == [option2]

    # we've only got chemistry, and no max
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    None, None) == [option1, option2]

def test_chem_bio_chem_steps_filtering():
    pathway = Pathway([Reaction('CCCO', ['CCC=O'], rxn_domain='chemistry'),
                       Reaction('CCC=O', ['CCC(=O)O'], rxn_domain='biocatalysis'),
                       Reaction('CCC(=O)O', ['CCC(=O)NC'], rxn_domain='chemistry')])

    option1 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'chemistry', 0.0, lambda x: x)
    option2 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'biocatalysis', 0.0, lambda x: x)

    options = [option1, option2]

    # no filtering because only 2 chemistry steps
    assert filter_by_max_chemistry_nodes(options, pathway, 4) == [option1, option2]

    # filter because max 1 step
    assert filter_by_max_chemistry_nodes(options, pathway, 1) == [option2]

    # we're at beginning, so can only add chemistry
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    2, 2) == [option1]

    # we're at beginning, but have reached max chem
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    1, 1) == []


def test_no_steps_filtering():
    pathway = Pathway([], target_smi='CCC(=O)NC')

    option1 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'chemistry', 0.0, lambda x: x)
    option2 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'biocatalysis', 0.0, lambda x: x)

    options = [option1, option2]

    # no filtering because no rxns
    assert filter_by_max_chemistry_nodes(options, pathway, 4) == [option1, option2]

    # no filtering because no rxns
    assert apply_beginning_and_end_chemistry_filter(options,
                                                   pathway,
                                                   1, 1) == [option1, option2]

def test_bio_bio_filtering():
    pathway = Pathway([Reaction('CCCO', ['CCC=O'], rxn_domain='biocatalysis'),
                       Reaction('CCC=O', ['CCC(=O)NC'], rxn_domain='biocatalysis')])

    option1 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'chemistry', 0.0, lambda x: x)
    option2 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'biocatalysis', 0.0, lambda x: x)

    options = [option1, option2]

    # no filtering because only biocatalysis
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    1, 1) == [option1, option2]


def test_chem_filtering():
    pathway = Pathway([Reaction('CCCO', ['CCC(=O)NC'], rxn_domain='chemistry')])

    option1 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'chemistry', 0.0, lambda x: x)
    option2 = ReactionOption('CCC(=O)NC', 'test', [], {}, 'test', 'biocatalysis', 0.0, lambda x: x)

    mcts_config = MCTS_Config()
    mcts_config.chemistry_only_at_beginning_or_end = True
    mcts_config.max_chemistry_at_end = None
    mcts_config.max_chemistry_at_beginning = None
    mcts_config.max_chemistry_nodes = 4

    options = [option1, option2]

    # no filtering because not at max yet
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    3, 3) == [option1, option2]

    # only biocatalysis possible because at max
    assert apply_beginning_and_end_chemistry_filter(options,
                                                    pathway,
                                                    1, 1) == [option2]

