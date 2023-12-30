

class SourceMol_Config():

    def __init__(self):
        self.target_always_not_buyable = True
        self.source_mol_mode = 'building_blocks'
        self.source_mol_commercial_vendors = ['sigma', 'apollo', 'flurochem', 'alfa', 'lifechem', 'molport']
        self.max_price_per_gram = 250
        self.block_if_price_over_max = True
        self.source_mol_metabolism_vendors = ['ecmdb']
        self.source_mols_can_be_chiral = True

    def set_to_all_commercial(self):
        self.source_mol_mode = 'building_blocks'
        self.source_mols_can_be_chiral = True
        self.source_mol_commercial_vendors = ['zinc', 'mcule', 'molport', 'sigma', 'apollo', 'flurochem', 'alfa', 'lifechem']

    def set_to_all_metabolites(self):
        self.source_mol_mode = 'metabolites'
        self.source_mols_can_be_chiral = True
        self.source_mol_metabolism_vendors = ['ecmdb']

    def get_mode_and_vendors(self):
        mode = self.source_mol_mode
        if mode not in ['building_blocks', 'metabolites']:
            Exception(f'Mode specified ({mode}) is invalid')
            return 0, {}

        vendors = None
        if mode == 'building_blocks':
            vendors = self.source_mol_commercial_vendors
        elif mode == 'metabolites':
            vendors = self.source_mol_metabolism_vendors

        return mode, vendors

    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

def source_mol_config_from_dict(config_dict: dict) -> SourceMol_Config:
    return SourceMol_Config().update_from_dict(config_dict)