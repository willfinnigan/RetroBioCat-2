from rdkit.Chem import rdChemReactions

from rbc2.template_application.apply_template.components.multi_step_applicator import Multi_step_applicator
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.components.rdkit_applicator import RdKit_Applicator
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction


class Test_multi_step_applicator():

    def test_multi_step_with_rdkit(self):
        multi_applicator = Multi_step_applicator()
        rdkit_applicator = RdKit_Applicator()
        rxn1 = rdChemReactions.ReactionFromSmarts("[#6:1]=[O:2]>>[#6:1]-[OH:2]")
        rxn2 = rdChemReactions.ReactionFromSmarts("[#6:1]-[O:2]>>[#6:1](=O)-[OH:2]")
        multi_rxns = {'test_1': [[[rxn1], [rxn2]]]}
        result = multi_applicator.apply_multi_step_rules(rdkit_applicator, 'C=O', multi_rxns)
        assert result['test_1'] == [['O=CO']]

    def test_multi_step_with_rdchiral(self):
        multi_applicator = Multi_step_applicator()
        rdchiral_applicator = RdChiral_Applicator()
        rxn1 = rdchiralReaction("[#6:1]=[O:2]>>[#6:1]-[OH:2]")
        rxn2 = rdchiralReaction("[#6:1]-[O:2]>>[#6:1](=O)-[OH:2]")
        multi_rxns = {'test_1': [[[rxn1], [rxn2]]]}
        result = multi_applicator.apply_multi_step_rules(rdchiral_applicator, 'C=O', multi_rxns)
        assert result['test_1'] == [['O=CO']]



