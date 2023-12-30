

class Expansion_Config():

    def __init__(self):

        # rule application
        self.allow_chiral_symmetry = False
        self.check_chiral_products = True
        self.combine_enantiomers = True
        self.allow_cyclic_reaction_outcomes = False
        self.clean_brackets = True

        # reaction parsing
        self.allow_backwards = False
        self.allow_duplicates = False
        self.duplicates_require_same_expander = True
        self.duplicates_require_same_domain = False
        self.duplicates_require_same_name = False
        self.merge_duplicate_metadata = True

        # expanders general
        self.max_reactions = None  # max reactions (not options)

        # reaction filtering and blocking
        self.use_max_mw_for_enzymes = False
        self.max_mw_to_use_enzymes = 300



    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

def expansion_config_from_dict(config_dict: dict) -> Expansion_Config:
    return Expansion_Config().update_from_dict(config_dict)



if __name__ == '__main__':
    config = Expansion_Config()
    dict_attr = config.to_dict()
    new_config = Expansion_Config().update_from_dict(dict_attr)
