from collections import defaultdict
from typing import Optional, List

from rdkit.Chem import rdChemReactions

from rbc2.template_application.apply_template.components.multi_step_applicator import Multi_step_applicator
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.components.rdkit_applicator import RdKit_Applicator
from rbc2.template_application.apply_template.components.template_flag_checks import is_reaction_cyclic, \
    is_intramolecular, is_dimer, does_reaction_keep_same_number_of_rings
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config


class RuleApplicator():

    def __init__(self, config=None):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.logger = add_logger('RuleApplicator', level=self.config.rule_applicator_log_level)

        self.rdkit_applicator = RdKit_Applicator(config=self.config, log_level=self.config.rule_applicator_log_level)
        self.rdchiral_applicator = RdChiral_Applicator(config=self.config, log_level=self.config.rule_applicator_log_level)
        self.multi_step_applicator = Multi_step_applicator(config=self.config, log_level=self.config.rule_applicator_log_level)
        self.logger = add_logger('RuleApplicator', level=self.config.rule_applicator_log_level)
        self.rule_applications = 0

    def apply(self, smi: str, rxns: dict[str: List],
              multistep_rxns=None, template_flags: Optional[dict] = None,
              return_atom_mapping=False):

        self.rule_applications += sum([len(v) for v in rxns.values()])

        if self.config.use_rdchiral == True:
            applicator = self.rdchiral_applicator
        else:
            applicator = self.rdkit_applicator

        if return_atom_mapping and self.config.use_rdchiral:
            product_dict, atom_mappings = applicator.apply_rules(smi, rxns,
                                                                  return_atom_mapping=True)
        else:
            product_dict = applicator.apply_rules(smi, rxns)
            atom_mappings = {}
            if return_atom_mapping and not self.config.use_rdchiral:
                self.logger.warning("return_atom_mapping is only supported with rdchiral; returning empty atom_mappings")

        if multistep_rxns is not None:
            if return_atom_mapping:
                multi_product_dict, multi_mapping_dict = self.multi_step_applicator.apply_multi_step_rules(
                    applicator, smi, multistep_rxns, return_atom_mapping=True)
                atom_mappings.update(multi_mapping_dict)
            else:
                multi_product_dict = self.multi_step_applicator.apply_multi_step_rules(
                    applicator, smi, multistep_rxns)
            product_dict.update(multi_product_dict)

        if return_atom_mapping:
            product_dict, atom_mappings = self._apply_template_flag_checks_with_mappings(
                smi, product_dict, template_flags, atom_mappings)
            return product_dict, atom_mappings

        product_dict = self._apply_template_flag_checks(smi, product_dict, template_flags)
        return product_dict

    # def _get_correct_format_rxns(self, rxns: List):
    #     """
    #     look at the first rxn - if its correct assume all are. If not, convert all to correct format
    #     rxns could be smarts, could be rdchiral reactions, or rdkit reactions."""
    #
    #     first_rxn = rxns[0]
    #     if isinstance(first_rxn, str):
    #         if self.config.use_rdchiral == True:
    #             self.logger.debug('Converting smarts to rdchiral')
    #             rxns = self.smarts_to_rdchiral(rxns)
    #         else:
    #             self.logger.debug('Converting smarts to rdkit')
    #             rxns = self.smarts_to_rdkit(rxns)
    #
    #     elif isinstance(first_rxn, rdchiralReaction):
    #         if self.config.use_rdchiral == False:
    #             self.logger.debug('Converting rdchiral to rdkit')
    #             rxns = [rxn.rxn for rxn in rxns]  # extracts the rdkit reaction from rdchiral
    #
    #     elif isinstance(first_rxn, rdChemReactions.ChemicalReaction):
    #         if self.config.use_rdchiral == True:
    #             self.logger.debug('Converting rdkit to rdchiral')
    #             smarts = [rxn.ToSmarts() for rxn in rxns]
    #             rxns = self.smarts_to_rdchiral(smarts)
    #
    #     return rxns



    def _apply_template_flag_checks_with_mappings(self, target_smi: str, product_dict: dict[str: list],
                                                      template_flags: Optional[dict], atom_mappings: dict):
        checked_product_dict = defaultdict(list)
        checked_mappings = defaultdict(list)

        if template_flags is None:
            template_flags = {}

        for name, outcome_list in product_dict.items():
            flags = template_flags.get(name, None)
            mapping_list = atom_mappings.get(name, [])
            for i, rxn_outcome in enumerate(outcome_list):
                if self._does_rxn_outcome_pass_checks(target_smi, rxn_outcome, flags) is True:
                    checked_product_dict[name].append(rxn_outcome)
                    if i < len(mapping_list):
                        checked_mappings[name].append(mapping_list[i])

        return checked_product_dict, dict(checked_mappings)

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

    applier = RuleApplicator()

    rxn1 = '[#6X4;z1:2]-[#7X3;z0:3]>>[#6X3;z1:2]=[#7X2;z0:3]'
    rxn2 = '[#6:1]-[#6@@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxn3 = '[#6:1]-[#6@H;X4:3](-[#6:2])-[#7X3;z0:4]>>[#6:1]-[#6H0;X3:3](-[#6:2])=[#7X2;z0:4]'
    rxns = {'Imine reduction': applier.smarts_to_rdchiral([rxn1, rxn2, rxn3])}

    products = applier.apply_rdchiral('[C@H]1(C2=CC=CC=C2)NCCCC1', rxns)

    print(products)
