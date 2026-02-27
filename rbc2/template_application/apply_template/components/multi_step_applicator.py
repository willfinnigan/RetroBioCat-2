from rbc2.configs.expansion_config import Expansion_Config
from rbc2.utils.add_logger import add_logger


class Multi_step_applicator():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.logger = add_logger('Multi_step_applicator', level=self.config.rule_applicator_log_level)

    def apply_multi_step_rules(self, applicator, smi, multi_step_rxns,
                                return_atom_mapping=False):
        rxn_product_dict = {}
        rxn_mapping_dict = {}
        for rxn_name, rxn_data in multi_step_rxns.items():
            self.logger.debug(f'Applying multi-step rule: {rxn_name}')

            if isinstance(rxn_data, tuple):
                grouped_list_steps, step_names_per_group = rxn_data
            else:
                grouped_list_steps = rxn_data
                step_names_per_group = None

            for group_idx, list_steps in enumerate(grouped_list_steps):
                if return_atom_mapping:
                    step_names = step_names_per_group[group_idx] if step_names_per_group else []
                    step_products, chains = self._apply_steps_with_mappings(
                        applicator, [smi], list_steps, step_names)
                else:
                    step_products = self._apply_steps(applicator, [smi], list_steps)
                    chains = []

                self.logger.debug(f'Products: {step_products}')

                if len(step_products) != 0:
                    rxn_product_dict.setdefault(rxn_name, []).extend(step_products)
                    if return_atom_mapping:
                        rxn_mapping_dict.setdefault(rxn_name, []).extend(chains)

        if return_atom_mapping:
            return rxn_product_dict, rxn_mapping_dict
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

    def _apply_steps_with_mappings(self, applicator, step_reactants, list_steps, step_names):
        """Apply sequential steps, building a mapping chain for each final product.

        When an intermediate step produces multiple substrates (e.g. ['A', 'B']),
        both substrates share the same chain prefix. This is a known approximation —
        in practice, multi-step reactions in the dataset produce single intermediates.
        """
        current = [(smi, []) for smi in step_reactants]

        for step_idx, step in enumerate(list_steps):
            name = step_names[step_idx] if step_idx < len(step_names) else ''
            next_current = []
            final_results = []

            for smi, chain in current:
                reactant = applicator.get_reactant(smi)
                products, mappings = applicator.apply_reactions(
                    reactant, step, return_atom_mapping=True)

                for product_set, mapping in zip(products, mappings):
                    step_dict = {
                        'rxn_name': name,
                        'target_smi': smi,
                        'substrate_smis': product_set,
                        'mapping': mapping
                    }
                    new_chain = chain + [step_dict]

                    if step_idx == len(list_steps) - 1:
                        final_results.append((product_set, new_chain))
                    else:
                        for sub_smi in product_set:
                            next_current.append((sub_smi, new_chain))

            if step_idx < len(list_steps) - 1:
                current = next_current
            else:
                return ([ps for ps, _ in final_results],
                        [ch for _, ch in final_results])

        return [], []
