from typing import List, Optional

from rdkit import Chem

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.template_application.apply_template.components.result_parser import Result_Parser
from rbc2.template_application.apply_template.rdchiral.clean import combine_enantiomers_into_racemic
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

    def apply_rules(self, smi: str, rxns: dict[str: List[rdchiralReaction]],
                    return_atom_mapping=False):
        reactant = self.get_reactant(smi)
        if reactant is None:
            return ({}, {}) if return_atom_mapping else {}

        rxn_products = {}
        rxn_mappings = {}
        for rxn_name, rxn_list in rxns.items():
            self.logger.debug(f'-- Applying {rxn_name} --')
            if return_atom_mapping:
                products, mappings = self.apply_reactions(reactant, rxn_list,
                                                          return_atom_mapping=True)
            else:
                products = self.apply_reactions(reactant, rxn_list)
                mappings = []
            if len(products) != 0:
                rxn_products[rxn_name] = products
                if return_atom_mapping:
                    rxn_mappings[rxn_name] = mappings

        if return_atom_mapping:
            return rxn_products, rxn_mappings
        return rxn_products

    def apply_reactions(self, reactant: rdchiralReactants, list_rdchiral_rxns: List[rdchiralReaction],
                        return_atom_mapping=False):
        if reactant is None:
            return ([], []) if return_atom_mapping else []

        all_products = []
        all_mappings = []
        for rdchiral_rxn in list_rdchiral_rxns:
            if return_atom_mapping:
                products, atom_mappings = self._apply_reaction(reactant, rdchiral_rxn,
                                                                return_atom_mapping=True)
                for p in products:
                    all_mappings.append(atom_mappings.get(p, {}))
            else:
                products = self._apply_reaction(reactant, rdchiral_rxn)
            all_products.extend(products)

        if len(all_products) == 0:
            return ([], []) if return_atom_mapping else []

        if return_atom_mapping:
            parsed, parsed_mappings = self._parse_with_mappings(all_products, all_mappings)
            return parsed, parsed_mappings
        else:
            parsed = self.parser.parse(all_products)
            return parsed

    def _parse_with_mappings(self, products, mappings):
        """Parse products while carrying atom mappings in parallel."""

        # combine enantiomers
        if self.config.combine_enantiomers:
            products, mappings = self._combine_enantiomers_with_mappings(products, mappings)

        # split
        split_products = self.parser.split_any_joined_molecules(products)

        # clean brackets
        if self.config.clean_brackets:
            split_products = self.parser.clean_brackets(split_products)

        # rdkit canonicalize + drop nones
        if self.config.force_rdkit_smis:
            split_products = self.parser.ensure_rdkit_products(split_products)
            kept_products, kept_mappings = [], []
            for sp, m in zip(split_products, mappings):
                if None not in sp:
                    kept_products.append(sp)
                    kept_mappings.append(m)
            split_products, mappings = kept_products, kept_mappings

        # dedup (keep first mapping for each unique substrate set)
        combined_products, combined_mappings = [], []
        for sp, m in zip(split_products, mappings):
            sp_set = list(set(sp))
            if sp_set not in combined_products:
                combined_products.append(sp_set)
                combined_mappings.append(m)

        return combined_products, combined_mappings

    def _combine_enantiomers_with_mappings(self, products, mappings):
        """Combine enantiomeric products while preserving parallel mappings.

        Runs the standard combine_enantiomers_into_racemic, then matches
        each resulting SMILES back to the first original that produced it
        by comparing stereo-stripped canonical forms.

        When two enantiomers merge into one achiral product, the first
        enantiomer's mapping is kept. This is safe because enantiomers from
        the same template match always map to the same target atoms — only
        the #3/#4 flanking-group assignment differs, which is meaningless
        once the stereocenter is removed. Products from different match
        positions always have structurally different achiral SMILES, so the
        stereo-stripped comparison is unambiguous across positions.
        """
        original_set = set(products)
        combined_set = combine_enantiomers_into_racemic(set(products))

        if original_set == combined_set:
            return products, mappings

        removed = original_set - combined_set
        added = combined_set - original_set

        # For each new (achiral) SMILES, find which removed original it came from
        # by comparing stereo-stripped canonical forms
        added_to_mapping = {}
        used_originals = set()
        for new_smi in added:
            new_mol = Chem.MolFromSmiles(new_smi)
            if new_mol:
                Chem.RemoveStereochemistry(new_mol)
                new_stripped = Chem.MolToSmiles(new_mol, True)
            else:
                new_stripped = new_smi

            for smi, m in zip(products, mappings):
                if smi in removed and smi not in used_originals:
                    old_mol = Chem.MolFromSmiles(smi)
                    if old_mol:
                        Chem.RemoveStereochemistry(old_mol)
                        old_stripped = Chem.MolToSmiles(old_mol, True)
                    else:
                        old_stripped = smi

                    if old_stripped == new_stripped:
                        added_to_mapping[new_smi] = m
                        used_originals.add(smi)
                        break

        # Build new parallel lists: survived originals (in order), then added
        result_products = []
        result_mappings = []
        seen = set()

        for smi, m in zip(products, mappings):
            if smi in combined_set and smi not in seen:
                result_products.append(smi)
                result_mappings.append(m)
                seen.add(smi)

        for new_smi in added:
            if new_smi not in seen:
                result_products.append(new_smi)
                result_mappings.append(added_to_mapping.get(new_smi, {}))
                seen.add(new_smi)

        return result_products, result_mappings

    def _apply_reaction(self, reactant: rdchiralReactants, rxn: rdchiralReaction,
                        return_atom_mapping=False):
        if reactant is None:
            return ([], {}) if return_atom_mapping else []

        try:
            if return_atom_mapping:
                result, atom_mappings = rdchiralRun(rxn, reactant,
                                                    combine_enantiomers=False,
                                                    check_chiral_products=self.config.check_chiral_products,
                                                    allow_chiral_symmetry=self.config.allow_chiral_symmetry,
                                                    return_atom_mapping=True)
            else:
                result = rdchiralRun(rxn, reactant,
                                     combine_enantiomers=False,
                                     check_chiral_products=self.config.check_chiral_products,
                                     allow_chiral_symmetry=self.config.allow_chiral_symmetry)
            if len(result) != 0:
                self.logger.debug(f'Applying rule {rxn.reaction_smarts} = {result}')
            if return_atom_mapping:
                return result, atom_mappings
            return result

        except Exception as e:
            self.logger.warning(f"Rdchiral rule applicator: Error running reactants for: {str(rxn.reaction_smarts)}")
            self.logger.warning(e)
            return ([], {}) if return_atom_mapping else []

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
