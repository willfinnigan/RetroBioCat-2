from rbc2.expansion.default_expander_interface import Expander
from rbc2.expansion.multi_expander import MultiExpander
from rbc2.expansion.expander_repository import get_expanders
from rbc2.data_model.reaction_option import ReactionOption


def test_mcts_expander_initialises_with_expanders():
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders)
    first_expander = list(multi_expander.expanders.values())[0]
    assert isinstance(first_expander, Expander)

def test_multi_expander_can_get_options_for_a_pathway_ordered_by_score():
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders)
    options = multi_expander.get_options(['C(N)C(N)CC=O'])
    multi_expander.expander_config.option_combination_method = 'order_by_score'

    for opt in options:
        print(f"{opt.rxn_type} - {opt.score}")

    assert isinstance(options[0], ReactionOption)
    assert options[0].score > options[1].score
    assert options[1].score > options[2].score
    assert options[2].score > options[3].score
    assert options[3].score > options[4].score
    assert options[4].score > options[5].score

def test_multi_expander_can_get_options_for_a_pathway_interleaved():
    expanders = get_expanders(('retrobiocat', 'aizynthfinder'))
    multi_expander = MultiExpander(expanders)
    options = multi_expander.get_options(['C(N)C(N)CC=O'], combination_by='interleave')

    assert options[0].rxn_type != options[1].rxn_type
    assert options[1].rxn_type != options[2].rxn_type
    assert options[2].rxn_type != options[3].rxn_type
    assert options[3].rxn_type != options[4].rxn_type
    assert options[4].rxn_type != options[5].rxn_type