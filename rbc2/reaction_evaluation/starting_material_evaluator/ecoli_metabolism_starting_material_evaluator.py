from functools import lru_cache
from typing import Optional, List

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_source_mols import does_source_mols_db_exist, download_source_mols_db
from rbc2.reaction_evaluation.starting_material_evaluator.sqlite_source_mol.connect_sqlitedb import SQLite_Database
from rbc2.reaction_evaluation.starting_material_evaluator.sqlite_source_mol.query_sqlitedb import DB_Query_SQLite
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluator

data_folder = f'{path_to_data_folder}/buyability'


class EcoliSME(StartingMaterialEvaluator):
    vendor_urls = {'ecmdb': 'https://ecmdb.ca/compounds/[[ID]]'}

    def __init__(self,
                 target_always_not_buyable: bool = True,
                 custom_smiles: List = None,
                 blocked_smiles: List = None,
                 source_mols_can_be_chiral: bool = True):

        self.target_always_not_buyable = target_always_not_buyable
        self.source_mol_metabolism_vendors = ['ecmdb']
        self.source_mols_can_be_chiral = source_mols_can_be_chiral

        if does_source_mols_db_exist() == False:
            download_source_mols_db()

        db_path = data_folder + '/source_mols.db'
        self.database = SQLite_Database(db_path)
        self.query = DB_Query_SQLite(self.database)
        self.cache_column_names = {}
        self.cache_vendor_names = {}

        self.custom_smiles = custom_smiles or []
        self.blocked_smiles = blocked_smiles or []

    @lru_cache(maxsize=10000)
    def eval(self, smi):
        if smi in self.blocked_smiles:
            return False, {}
        if smi in self.custom_smiles:
            return True, {}

        mode = 'metabolites'
        vendors = self.source_mol_metabolism_vendors

        if self.is_mol_chiral(smi) and self.source_mols_can_be_chiral is False:
            return False, {}

        result = self.query.smiles_lookup(smi, mode, vendors=vendors)
        if result is None:
            return False, {}

        info = self._process_info(result, mode)

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



if __name__ == '__main__':
    sme = EColiSME()
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