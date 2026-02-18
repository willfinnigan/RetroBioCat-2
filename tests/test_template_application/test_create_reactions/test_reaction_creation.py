
from rbc2.data_model.reaction import Reaction
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions


def test_reaction_is_created_from_rule_application_output():
    output = {'test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', output, rxn_type='test', rxn_domain='not real')
    assert len(reactions) == 1
    assert reactions[0].substrates == ['CCO']
    assert reactions[0].product == 'CC=O'

def test_fwd_reaction_with_multiple_products():
    output = {'test_rxn': [['CCO', 'CC']]}
    reactions = create_reactions('CC=O.C', output, rxn_type='test', rxn_domain='not real', fwd_rxn=True)
    assert len(reactions) == 1
    assert reactions[0].product == 'CCO.CC'
    assert reactions[0].substrates == ['CC=O', 'C']

def test_reactions_without_init_uuid_have_different_uuid():
    rxn1 = Reaction(product='CC=O', substrates=['CCO'], rxn_type='test', rxn_domain='not real')
    rxn2 = Reaction(product='CC=O', substrates=['CCO'], rxn_type='test', rxn_domain='not real')
    assert rxn1.unique_id != rxn2.unique_id
