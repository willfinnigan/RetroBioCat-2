from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import reaction_from_dict
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions
from rbc2.template_application.create_reactions_from_output.find_duplicates import remove_duplicates
from rbc2.template_application.create_reactions_from_output.process_reactions import process_reactions


def make_test_network():
    network = Network()
    rule_applicator_output = {'test_rxn': [['CCO']]}
    reactions = create_reactions(target_smi='CC=O',
                                 outcomes_dict=rule_applicator_output,
                                 metadata=None,
                                 rxn_type='test',
                                 rxn_domain='not real')
    network.add_reaction(reactions[0])
    return network


def test_duplicate_reactions_are_found_if_identical_smiles():

    network = make_test_network()

    example_reaction_outcome = {'matching_test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome, rxn_type='different_rxn_type', rxn_domain='different_domain')
    reactions = remove_duplicates(reactions,
                                  network,
                                  require_same_expander=False,
                                  require_same_domain=False,
                                  require_matching_name=False,
                                  update_match_metadata=True)
    assert len(reactions) == 0


def test_duplicate_reactions_not_found_if_domain_is_different():

    network = make_test_network()

    example_reaction_outcome = {'matching_test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome, rxn_type='different_rxn_type', rxn_domain='different_domain')
    reactions = remove_duplicates(reactions,
                                  network,
                                  require_same_expander=False,
                                  require_same_domain=True,
                                  require_matching_name=False,
                                  update_match_metadata=True)
    assert len(reactions) == 1

def test_duplicate_reactions_not_found_if_rxn_type_is_different():

    network = make_test_network()

    example_reaction_outcome = {'matching_test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome, rxn_type='different_rxn_type', rxn_domain='different_domain')
    reactions = remove_duplicates(reactions,
                                  network,
                                  require_same_expander=True,
                                  require_same_domain=False,
                                  require_matching_name=False,
                                  update_match_metadata=True)
    assert len(reactions) == 1

def test_duplicate_reactions_are_found_if_identical_smiles_from_same_expander():

    network = make_test_network()

    example_reaction_outcome = {'test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome, rxn_type='test', rxn_domain='different_domain')
    assert len(reactions) != 0
    reactions = remove_duplicates(reactions,
                                  network,
                                  require_same_expander=True,
                                  require_same_domain=False,
                                  require_matching_name=True,
                                  update_match_metadata=False)
    assert len(reactions) == 0

def test_duplicate_reactions_not_found_if_rxn_name_is_different():

    network = make_test_network()

    example_reaction_outcome = {'matching_test_rxn': [['CCO']]}
    reactions = create_reactions('CC=O', example_reaction_outcome, rxn_type='different_rxn_type', rxn_domain='different_domain')
    reactions = remove_duplicates(reactions,
                                  network,
                                  require_same_expander=False,
                                  require_same_domain=False,
                                  require_matching_name=True,
                                  update_match_metadata=False)
    assert len(reactions) == 1

def test_duplicate_detection_of_process_reactions():
    reaction_dicts = [{'product': '[C@H]1(C2=CC=CC=C2)NCCCC1', 'substrates': ['C1=N[C@H](c2ccccc2)CCC1'],
                       'unique_id': 'Imine reduction_id_ca4b38f8-4727-46ff-ac82-d364eea11cb9',
                       'name': 'Imine reduction', 'rxn_type': 'retrobiocat', 'rxn_domain': 'biocatalysis'},
                      {'product': '[C@H]1(C2=CC=CC=C2)NCCCC1', 'substrates': ['c1ccc(C2=NCCCC2)cc1'],
                       'unique_id': 'Imine reduction_id_56210a51-3943-4ab7-9e9d-683d393220bd',
                       'name': 'Imine reduction', 'rxn_type': 'retrobiocat', 'rxn_domain': 'biocatalysis'}]
    local_reactions = [reaction_from_dict(reaction_dict) for reaction_dict in reaction_dicts]
    network = Network(reactions=local_reactions)

    test_reaction_dict = {'product': '[C@H]1(C2=CC=CC=C2)NCCCC1', 'substrates': ['C1=N[C@H](c2ccccc2)CCC1'],
                         'unique_id': 'different_unique_id',
                         'name': 'Imine reduction', 'rxn_type': 'retrobiocat', 'rxn_domain': 'biocatalysis'}
    reactions = [reaction_from_dict(test_reaction_dict)]

    reactions = process_reactions(reactions, network, Expansion_Config())

    assert len(reactions) == 0



