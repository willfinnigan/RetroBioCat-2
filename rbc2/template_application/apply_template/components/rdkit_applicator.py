from itertools import permutations

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from rbc2.template_application.apply_template.components.result_parser import Result_Parser
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config


class RdKit_Applicator():

    def __init__(self, config=None, log_level='WARNING'):

        self.config = config
        if self.config is None:
            self.config = Expansion_Config()
        self.parser = Result_Parser(config=self.config, log_level=log_level)
        self.logger = add_logger('RDKitApplicator', level=log_level)

    def apply_rules(self, smi: str, rxns: dict):
        reactant = self.get_reactant(smi)
        if reactant is None:
            return {}

        rxn_products = {}
        for rxn_name, rxn_list in rxns.items():
            products = self.apply_reactions(reactant, rxn_list)
            if len(products) != 0:
                rxn_products[rxn_name] = products

        return rxn_products

    def _apply_reaction(self, reactant: list[Chem.Mol], rxn: rdChemReactions.ChemicalReaction):
        if reactant is None:
            return []

        try:
            if len(reactant) > 1:
                products = []
                for perm in permutations(reactant):
                    products.extend(rxn.RunReactants(tuple(perm)))
            else:
                products = rxn.RunReactants(tuple(reactant))
            products = self.rdkit_out_to_rdchiral_output(products)
            return products

        except Exception as e:
            self.logger.error(f"Rdkit rule applicator: Error running reactants")
            self.logger.error(e)
            return []

    def apply_reactions(self, reactant: list[Chem.Mol], list_rdkit_rxns):
        if reactant is None:
            return []

        all_products = []
        for rdkit_rxn in list_rdkit_rxns:
            all_products.extend(self._apply_reaction(reactant, rdkit_rxn))

        if len(all_products) == 0:
            return []

        return self.parser.parse(all_products)

    def rdkit_out_to_rdchiral_output(self, rdkit_reaction_products):
        outcomes = []
        for products in rdkit_reaction_products:
            smi = ""
            for p in products:
                s = Chem.MolToSmiles(p)
                if smi == "":
                    smi = s
                else:
                    smi += f".{s}"
            outcomes.append(smi)
        return outcomes

    def get_reactant(self, smile_string):
        try:
            smis = smile_string.split('.')
            return [Chem.MolFromSmiles(s) for s in smis]

        except Exception as e:
            self.logger.error(f"Error creating mol for smiles string: {smile_string}")
            self.logger.error(e)
            return None

    def get_rxns(self, list_smarts):
        list_rxns = []
        for smarts in list_smarts:
            try:
                rxn = rdChemReactions.ReactionFromSmarts(smarts)
                list_rxns.append(rxn)
            except Exception as e:
                self.logger.error(f"Error creating rxn for smarts string: {smarts}")
                self.logger.error(e)
                return None

        return list_rxns




if __name__ == '__main__':
    applicator = RdKit_Applicator()
    rxn = rdChemReactions.ReactionFromSmarts("[#6:1]=[O:2]>>[#6:1]-[OH:2]")
    reactant = Chem.MolFromSmiles('C=O')
    products = applicator._apply_reaction(reactant, rxn)
    for prods in products:
        for p in prods:
            print(Chem.MolToSmiles(p))

    products = applicator.apply_reactions(reactant, rxn)
    for prods in products:
        for p in prods:
            print(Chem.MolToSmiles(p))