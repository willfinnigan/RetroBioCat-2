from rbc2.configs.expansion_config import Expansion_Config

class Multi_step_applicator():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

    def apply_multi_step_rules(self, applicator, smi, multi_step_rxns):
        rxn_product_dict = {}
        for rxn_name, grouped_list_steps in multi_step_rxns.items():
            for list_steps in grouped_list_steps:
                step_products = self._apply_steps(applicator, [smi], list_steps)

                if len(step_products) != 0:
                    if rxn_name not in rxn_product_dict:
                        rxn_product_dict[rxn_name] = []
                    rxn_product_dict[rxn_name] += step_products

        return rxn_product_dict

    def _apply_steps(self, applicator, step_reactants, list_steps):
        step_products = []
        for step in list_steps:
            step_products = []
            for smi in step_reactants:
                reactant = applicator.get_reactant(smi)
                step_products.extend(applicator.apply_reactions(reactant, step))

            step_reactants = []
            for smi_list in step_products:
                step_reactants.extend(smi_list)
        return step_products
