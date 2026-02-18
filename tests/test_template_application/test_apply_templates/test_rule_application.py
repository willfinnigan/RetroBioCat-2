import pytest

from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.components.rdkit_applicator import RdKit_Applicator
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator


@pytest.mark.parametrize('applicator', [RdChiral_Applicator(), RdKit_Applicator()])
class Test_rule_applicators():

    def test_apply_basic_oxidation_gives_expected_product(self, applicator):
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        reactant = applicator.get_reactant('C=O')
        products = applicator.apply_reactions(reactant, rxns)

        assert products == [['CO']]

    def test_apply_multiple_reactions(self, applicator):
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]", "[#6:1](-[#6:2])=[O:3]>>[#6:1](-[#6:2])-[OH:3]"])
        reactant = applicator.get_reactant('C=O')
        products = applicator.apply_reactions(reactant, rxns)

        assert products == [['CO']]

    def test_apply_dict_of_rules(self, applicator):
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]", "[#6:1](-[#6:2])=[O:3]>>[#6:1](-[#6:2])-[OH:3]"])
        rxns_dict = {'test_1': rxns}
        smi = 'C=O'
        products = applicator.apply_rules(smi, rxns_dict)

        assert products == {'test_1': [['CO']]}

    def test_apply_dict_of_rules_empty_rxns_not_returned(self, applicator):
        rxns = {'test_1': applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"]),
                'test_2': applicator.get_rxns(["[Cl:1]>>[Cl:1]-C"])}
        smi = 'C=O'
        products = applicator.apply_rules(smi, rxns)

        assert 'test_2' not in products


def test_rule_applicator_counts():
    applicator = RuleApplicator()
    rxns = {'test_1': applicator.smarts_to_rdchiral(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"]),
            'test_2': applicator.smarts_to_rdchiral(["[Cl:1]>>[Cl:1]-C"])}
    smi = 'C=O'
    products = applicator.apply(smi, rxns)
    assert applicator.rule_applications == 2

def test_multi_reactant_rule():
    applicator = RdKit_Applicator(log_level='DEBUG')
    rxns = applicator.get_rxns(["[#6:1]=[O:2].[#6:3]-[N:4]>>[#6:1]-[N:4]-[#6:3]"])
    reactant = applicator.get_reactant('C=O.CN')
    products = applicator.apply_reactions(reactant, rxns)

def test_multi_reactant_rule_with_reversed_order():
    """Reactants given in opposite order to SMARTS template should still produce results via permutation."""
    applicator = RdKit_Applicator()
    # SMARTS expects ketone first, amine second
    rxns = applicator.get_rxns(["[#6:1]=[O:2].[#6:3]-[N:4]>>[#6:1]-[N:4]-[#6:3]"])
    # Input has amine first, ketone second
    reactant = applicator.get_reactant('CN.C=O')
    products = applicator.apply_reactions(reactant, rxns)
    assert len(products) > 0