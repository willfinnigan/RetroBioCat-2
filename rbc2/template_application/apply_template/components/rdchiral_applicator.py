from typing import List, Optional

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.template_application.apply_template.components.result_parser import Result_Parser
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rbc2.template_application.apply_template.rdchiral.main import rdchiralRun
from rbc2.utils.add_logger import add_logger


class RdChiral_Applicator():

    def __init__(self, config=None, log_level='ERROR'):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.parser = Result_Parser(config=self.config, log_level=log_level)
        self.logger = add_logger('RDChiralApplicator', level=log_level)

    def apply_rules(self, smi: str, rxns: dict[str: List[rdchiralReaction]]):
        reactant = self.get_reactant(smi)
        if reactant is None:
            return {}

        rxn_products = {}
        for rxn_name, rxn_list in rxns.items():
            self.logger.debug(f'-- Applying {rxn_name} --')
            products = self.apply_reactions(reactant, rxn_list)
            if len(products) != 0:
                rxn_products[rxn_name] = products

        return rxn_products

    def apply_reactions(self, reactant: rdchiralReactants, list_rdchiral_rxns: List[rdchiralReaction]):
        if reactant is None:
            return []
        all_products = []
        for rdchiral_rxn in list_rdchiral_rxns:
            all_products.extend(self._apply_reaction(reactant, rdchiral_rxn))

        if len(all_products) == 0:
            return []

        return self.parser.parse(all_products)

    def _apply_reaction(self, reactant: rdchiralReactants, rxn: rdchiralReaction):
        if reactant is None:
            return []

        try:
            result = rdchiralRun(rxn, reactant,
                                 combine_enantiomers=False,
                                 check_chiral_products=self.config.check_chiral_products,
                                 allow_chiral_symmetry=self.config.allow_chiral_symmetry)
            if len(result) != 0:
                self.logger.debug(f'Applying rule {rxn.reaction_smarts} = {result}')
            return result

        except Exception as e:
            self.logger.warning(f"Rdchiral rule applicator: Error running reactants for: {str(rxn.reaction_smarts)}")
            self.logger.warning(e)
            return []

    def get_reactant(self, smile_string: str) -> Optional[rdchiralReactants]:

        if '.' in smile_string:
            self.logger.warning("RdChiral doesn't take multiple reactants, return empty result")
            return None

        try:
            return rdchiralReactants(smile_string)

        except Exception as e:
            self.logger.error(f"Error creating rchiral reaction for smiles string: {smile_string}")
            self.logger.error(e)
            return None

    def get_rxns(self, list_smarts, remove_incorrect_atom_numbering=False) -> List[rdchiralReaction]:
        list_rxns = []
        for smarts in list_smarts:
            try:
                rxn = rdchiralReaction(smarts, remove_incorrect_atom_numbering=remove_incorrect_atom_numbering)
                list_rxns.append(rxn)
            except Exception as e:
                self.logger.error(f"Error creating rxn for smarts string: {smarts}")
                self.logger.error(e)

        return list_rxns

