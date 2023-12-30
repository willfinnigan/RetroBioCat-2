from collections import defaultdict
from typing import Optional, List

from rbc2.template_application.apply_template.components.multi_step_applicator import Multi_step_applicator
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.components.rdkit_applicator import RdKit_Applicator
from rbc2.template_application.apply_template.components.template_flag_checks import is_reaction_cyclic, \
    is_intramolecular, is_dimer, does_reaction_keep_same_number_of_rings
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config


class RuleApplicator():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.rdkit_applicator = RdKit_Applicator(config=self.config, log_level=log_level)
        self.rdchiral_applicator = RdChiral_Applicator(config=self.config, log_level=log_level)
        self.multi_step_applicator = Multi_step_applicator(config=self.config, log_level=log_level)
        self.logger = add_logger('RuleApplicator', level=log_level)
        self.rule_applications = 0

    def apply_rdkit(self, smi: str, rxns: dict[str: List[str]],
                    multistep_rxns=None, template_flags: Optional[dict] = None):

        self.rule_applications += sum([len(v) for v in rxns.values()])
        product_dict = self.rdkit_applicator.apply_rules(smi, rxns)

        if multistep_rxns is not None:
            multi_product_dict = self.multi_step_applicator.apply_multi_step_rules(self.rdkit_applicator, smi, multistep_rxns)
            product_dict.update(multi_product_dict)

        product_dict = self._apply_template_flag_checks(smi, product_dict, template_flags)

        return product_dict

    def apply_rdchiral(self, smi: str, rxns: dict[str: List[rdchiralReaction]],
                       multistep_rxns=None, template_flags: Optional[dict] = None):

        self.rule_applications += sum([len(v) for v in rxns.values()])
        product_dict = self.rdchiral_applicator.apply_rules(smi, rxns)

        if multistep_rxns is not None:
            multi_product_dict = self.multi_step_applicator.apply_multi_step_rules(self.rdchiral_applicator, smi,
                                                                                   multistep_rxns)
            product_dict.update(multi_product_dict)

        product_dict = self._apply_template_flag_checks(smi, product_dict, template_flags)

        return product_dict

    def _apply_template_flag_checks(self, target_smi: str, product_dict: dict[str: list], template_flags: Optional[dict]):
        checked_product_dict = defaultdict(list)

        if template_flags is None:
            template_flags = {}

        for name, outcome_list in product_dict.items():
            flags = template_flags.get(name, None)
            for rxn_outcome in outcome_list:
                if self._does_rxn_outcome_pass_checks(target_smi, rxn_outcome, flags) is True:
                    checked_product_dict[name].append(rxn_outcome)

        return checked_product_dict

    def _does_rxn_outcome_pass_checks(self, target_smi, rxn_outcome, template_flags):

        if is_reaction_cyclic(target_smi, rxn_outcome):
            return False

        if template_flags is None:
            return True

        if template_flags.get('intra_only', False) is True:
            if is_intramolecular(rxn_outcome) is False:
                return False

        if template_flags.get('dimer_only', False) is True:
            if is_dimer(rxn_outcome) is False:
                return False

        if template_flags.get('ring_break', False) is True:
            if does_reaction_keep_same_number_of_rings(target_smi, rxn_outcome):
                return False

        return True

    def smarts_to_rdchiral(self, smarts, remove_incorrect_atom_numbering=False) -> List[rdchiralReaction]:
        return self.rdchiral_applicator.get_rxns(smarts, remove_incorrect_atom_numbering=remove_incorrect_atom_numbering)

    def smarts_to_rdkit(self, smarts):
        return self.rdkit_applicator.get_rxns(smarts)

if __name__ == '__main__':

    applier = RuleApplicator(log_level='DEBUG')

    rxn1 = '[#6X4;z1:2]-[#7X3;z0:3]>>[#6X3;z1:2]=[#7X2;z0:3]'
    rxn2 = '[#6:1]-[#6@@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxn3 = '[#6:1]-[#6@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxns = {'Imine reduction': applier.smarts_to_rdchiral([rxn1, rxn2, rxn3])}

    products = applier.apply_rdchiral('[C@H]1(C2=CC=CC=C2)NCCCC1', rxns)

    print(products)
