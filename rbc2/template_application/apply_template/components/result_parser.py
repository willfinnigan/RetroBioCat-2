from rbc2.template_application.apply_template.components.bracket_cleaner import BracketCleaner
from rbc2.template_application.apply_template.rdchiral.clean import combine_enantiomers_into_racemic
from rbc2.utils.add_logger import add_logger
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.utils.rdkit_utils import rdkit_smile


class Result_Parser():

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = Expansion_Config()

        self.bracket_cleaner = BracketCleaner(log_level=log_level)
        self.logger = add_logger("Result_Parser", level=log_level)

    def parse(self, reaction_products):

        self.logger.debug(f"Parsing reaction products: {reaction_products}")

        if self.config.combine_enantiomers:
            reaction_products = self.combine_enantiomers(reaction_products)

        reaction_products = self.split_any_joined_molecules(reaction_products)

        if self.config.clean_brackets:
            reaction_products = self.clean_brackets(reaction_products)

        if self.config.force_rdkit_smis:
            reaction_products = self.ensure_rdkit_products(reaction_products)
            reaction_products = self.drop_any_nones(reaction_products)

        reaction_products = self.combine_products(reaction_products)

        return reaction_products

    def drop_any_nones(self, parsed_reaction_products):
        kept_products = []
        for list_smis in parsed_reaction_products:
            if None not in list_smis:
                kept_products.append(list_smis)
            else:
                self.logger.debug(f"Dropped smis with none - {list_smis}")
        return kept_products

    def ensure_rdkit_products(self, parsed_reaction_products):
        processed = []
        for products in parsed_reaction_products:
            processed.append(self._rdkit_products(products))
        return processed

    def combine_products(self, reaction_products):

        reaction_product_sets = []
        for products in reaction_products:
            reaction_product_sets.append(set(products))

        combined = []
        for products in reaction_product_sets:
            list_products = list(products)
            if list_products not in combined:
                combined.append(list_products)

        if reaction_products != combined:
            self.logger.debug(f"Combined reaction products. Original = {reaction_products}, New = {combined}")

        return combined

    def combine_enantiomers(self, reaction_products):
        reaction_products = combine_enantiomers_into_racemic(set(reaction_products))
        self.logger.debug(f"Combined enantiomers: {reaction_products}")
        return reaction_products

    def split_any_joined_molecules(self, reaction_products):
        processed_smis = []
        for smi in reaction_products:
            processed_smis.append(self._split_products(smi))
        self.logger.debug(f"Split up smis:  {processed_smis}")
        return processed_smis

    def clean_brackets(self, parsed_reaction_products):
        cleaned_products = []
        for list_smi in parsed_reaction_products:
            cleaned_list = [self.bracket_cleaner.clean_brackets(smi) for smi in list_smi]
            cleaned_products.append(cleaned_list)

        if cleaned_products != parsed_reaction_products:
            self.logger.debug(f'-- Brackets cleaned --')
            self.logger.debug(f"Original: {parsed_reaction_products}")
            self.logger.debug(f"Cleaned: {cleaned_products}")
        return cleaned_products

    def _split_products(self, smi):
        split_list = smi.split('.')
        return split_list

    def _rdkit_products(self, listSmi):
        new_list = []
        for smi in listSmi:
            new_list.append(rdkit_smile(smi))
        return new_list

if __name__ == '__main__':
    parser = Result_Parser(log_level='DEBUG')
    #products = ['CO.C', 'CO.C']
    #parsed_products = parser.parse(products)
    #print(parsed_products)

    products = ['[C]1=[C]C[C@@H](c2ccccc2)N=C1', 'C1=CC[C@@H](c2ccccc2)N=C1']
    parsed_products = parser.parse(products)
    print(parsed_products)

