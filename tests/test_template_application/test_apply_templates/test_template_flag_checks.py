from rbc2.template_application.apply_template.rule_applicator import RuleApplicator


def test_rule_applicator_for_basic_oxidation():
    applicator = RuleApplicator()
    smarts = ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"]
    rxns = {'test_rxn': applicator.smarts_to_rdchiral(smarts)}
    flags = {'test_rxn': {'intra_only': True}}
    result = applicator.apply_rdchiral("CC=O", rxns, template_flags=flags)

    assert result == {'test_rxn': [['CCO']]}

def test_cyclic_products_are_removed():
    applicator = RuleApplicator()
    smarts = ["[#6:1]=[O:2]>>[#6:1]=[O:2]"]
    rxns = {'test_rxn': applicator.smarts_to_rdchiral(smarts)}
    flags = {'test_rxn': {'intra_only': True}}
    result = applicator.apply_rdchiral("CC=O", rxns, template_flags=flags)
    assert result == {}

def test_non_intra_products_are_removed():
    applicator = RuleApplicator()
    smarts = ["[#6:1]=[O:2]>>[#6:1]-[OH:2].CC"]
    rxns = {'test_rxn': applicator.smarts_to_rdchiral(smarts)}
    flags = {'test_rxn': {'intra_only': True}}
    result = applicator.apply_rdchiral("CC=O", rxns, template_flags=flags)
    assert result == {}


def test_non_dimer_products_are_removed():
    applicator = RuleApplicator()
    smarts = ["[#6:1]=[O:2]>>[#6:1]-[OH:2].CC"]
    rxns = {'test_rxn': applicator.smarts_to_rdchiral(smarts)}
    flags = {'test_rxn': {'dimer_only': True}}
    result = applicator.apply_rdchiral("CC=O", rxns, template_flags=flags)
    assert result == {}

