from functools import lru_cache
from typing import Optional

from rbc2.configs.data_path import path_to_data_folder
from rbc2.reaction_evaluation.starting_material_evaluator.sqlite_source_mol.connect_sqlitedb import SQLite_Database
from rbc2.reaction_evaluation.starting_material_evaluator.sqlite_source_mol.query_sqlitedb import DB_Query_SQLite
from rbc2.configs.source_mol_config import SourceMol_Config
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluatorInterface

data_folder = f'{path_to_data_folder}/buyability'

class DefaultSQLStartingMaterialEvaluator(StartingMaterialEvaluatorInterface):
    vendor_urls = {'mcule': 'https://mcule.com/[[ID]]',
                   'sigma': 'https://www.sigmaaldrich.com/GB/en/search/[[ID]]?focus=products&page=1&perpage=30&sort=relevance&term=[[ID]]&type=product',
                   'lifechem': 'https://shop.lifechemicals.com/compound/[[ID]]',
                   'apollo': 'https://store.apolloscientific.co.uk/search?search=[[ID]]',
                   'alfa': 'https://www.alfa.com/en/catalog/[[ID]]',
                   'zinc': 'https://zinc.docking.org/substances/[[ID]]',
                   'flurochem': 'http://www.fluorochem.co.uk/Products/Product?code=[[ID]]',
                   'molport': 'https://www.molport.com/shop/molecule-link/[[ID]]',
                   'ecmdb': 'https://ecmdb.ca/compounds/[[ID]]'}

    available_modes = ['building_blocks', 'metabolites']

    def __init__(self, config: Optional[SourceMol_Config] = None, custom_smiles=None, blocked_smiles=None):

        db_path = data_folder + '/source_mols.db'
        self.database = SQLite_Database(db_path)
        self.query = DB_Query_SQLite(self.database)
        self.cache_column_names = {}
        self.cache_vendor_names = {}
        self.config = config
        if self.config is None:
            self.config = SourceMol_Config()

        self.custom_smiles = custom_smiles
        if self.custom_smiles is None:
            self.custom_smiles = []

        self.blocked_smiles = blocked_smiles
        if self.blocked_smiles is None:
            self.blocked_smiles = []

    @lru_cache(maxsize=10000)
    def eval(self, smi):
        if smi in self.blocked_smiles:
            return False, {}
        if smi in self.custom_smiles:
            return True, {}

        mode, vendors = self.config.get_mode_and_vendors()

        if self.is_mol_chiral(smi) and self.config.source_mols_can_be_chiral is False:
            return False, {}

        result = self.query.smiles_lookup(smi, mode, vendors=vendors)
        if result is None:
            return False, {}

        info = self._process_info(result, mode)

        if self._is_above_max_price_per_gram(info, vendors) == True:
            return False, info
        return True, info

    def is_mol_chiral(self, smi):
        if '@' in smi:
            return True
        return False

    @lru_cache(maxsize=10)
    def column_names(self, mode):
        if mode not in self.cache_column_names:
            self.cache_column_names[mode] = self.query.get_column_names(mode)
        return self.cache_column_names[mode]

    @lru_cache(maxsize=10)
    def vendor_names(self, mode):
        if mode not in self.cache_vendor_names:
            columns = self.column_names(mode)
            vendors = []
            for col in columns:
                if '_id' in col:
                    vendors.append(col.replace('_id', ''))
            self.cache_vendor_names[mode] = vendors
        return self.cache_vendor_names[mode]

    def _process_info(self, result, mode):
        columns = self.column_names(mode)
        vendors = self.vendor_names(mode)
        info = {k: v for k, v in zip(columns, result)}
        info = {k: v for k, v in info.items() if v is not None}

        vendor_info = {}
        for col, value in info.items():
            for vendor in vendors:
                if vendor in col:
                    if vendor not in vendor_info:
                        vendor_info[vendor] = {}
                    vendor_info[vendor][col.replace(f"{vendor}_", '')] = value
                    if ('_id' in col) and (vendor in self.vendor_urls):
                        url = self.vendor_urls[vendor].replace('[[ID]]', value)
                        vendor_info[vendor]['url'] = url
        return vendor_info

    def _is_above_max_price_per_gram(self, info, requested_vendors):
        """ Determines whether the price of the molecule is above the maximum price per gram, based on settings in the config """

        if requested_vendors is None:
            return False

        price_too_high = []  # will become list of booleans, one for each vendor, for whether the price is too high
        for vendor in requested_vendors:
            if vendor in info:
                if 'ppg' in info[vendor]:
                    if info[vendor]['ppg'] is None:
                        price_too_high.append(False)
                    elif float(info[vendor]['ppg']) > self.config.max_price_per_gram:
                        price_too_high.append(True)
                    else:
                        price_too_high.append(False)
                else:
                    price_too_high.append(False)


        # if there are only Trues in the list, return True
        if len(price_too_high) == sum(price_too_high):
            return True

        # if config.block_if_price_over_max is True, then return True if any of the prices are too high
        if self.config.block_if_price_over_max is True:
            if True in price_too_high:
                return True

        # otherwise return False
        return False



if __name__ == '__main__':
    sme = DefaultSQLStartingMaterialEvaluator()
    available, info = sme.eval('CCCC=O')
    print(info)
    print(available)


    """
    import time

    for i in range(20):
        sme = StartingMaterialEvaluator()
        t0 = time.time()
        available, info = sme.eval('CCC(=O)C(=O)O')
        t1 = time.time()
        print(f'Time for 1 = {round((t1-t0)*1000,3)} milliseconds')

    def time_test():
        t0 = time.time()
        sme = StartingMaterialEvaluator()
        for i in range(20000):
            available, info = StartingMaterialEvaluator().eval('CCC(=O)C(=O)O')
            #available, info = sme.eval('CCC(=O)C(=O)O')
        t1 = time.time()
        time_taken = t1 - t0
        print(f"Total time = {time_taken}")
        print(f"Time per query = {time_taken / 10000}")

    time_test()
    """