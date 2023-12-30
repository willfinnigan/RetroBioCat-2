from rbc2.utils.add_logger import add_logger
import itertools

class BracketCleaner():
    """ This module cleans up unwanted brackets which makes matching smiles strings problematic """

    def __init__(self, log_level='WARNING'):
        self.bracket_required_chars = ['@', '+', '-', 'H', 'h']
        self.bracker_must_contain_for_cleaning = ['C', 'N', 'O']
        self.bracket_remove = ['[', ']']
        self.logger = add_logger("Bracket cleaner", level=log_level)

    def clean_brackets(self, smi):
        """
        If smi contains '[', check if brackets are needed, if not then remove.
        """

        if '[' in smi:
            bracket_locations = self._get_bracket_locations(smi)
            self.logger.debug(f"Bracket locations = {bracket_locations}")

            keep_brackets = self._get_locations_required(smi, bracket_locations)
            keep_locations = [loc for loc, keep in zip(bracket_locations, keep_brackets) if not keep]

            smi = self._remove_brackets(smi, keep_locations)

        return smi

    def _get_bracket_locations(self, smi):
        fwd_brackets = [i for i, l in enumerate(smi) if l == '[']
        rev_brackets = [i for i, l in enumerate(smi) if l == ']']
        locations = list(zip(fwd_brackets, rev_brackets))
        return locations

    def _get_locations_required(self, smi, locations):
        brackets_to_remove = []
        for loc in locations:
            bracket = smi[loc[0]+1:loc[1]]
            keep_bracket = False

            if self._does_bracket_contain_required_char(bracket) or not self._does_bracket_contain_char_to_warrent_cleaning(bracket):
                keep_bracket = True
                brackets_to_remove.append(True)
            brackets_to_remove.append(keep_bracket)
            self.logger.debug(f"Bracket contents at {loc} = {bracket} - (keep={keep_bracket})")

        return brackets_to_remove

    def _does_bracket_contain_required_char(self, bracket):
        for l in bracket:
            if l in self.bracket_required_chars:
                return True
        return False

    def _does_bracket_contain_char_to_warrent_cleaning(self, bracket):
        for l in self.bracker_must_contain_for_cleaning:
            if l in bracket:
                return True
        self.logger.debug(f"Bracket {bracket} does not contain a char required to warrent cleaning")
        return False

    def _remove_brackets(self, smi, locations):
        all_locations = list(itertools.chain(*locations))
        new_smi = ""
        for i, l in enumerate(smi):
            if i not in all_locations:
                new_smi += l
        self.logger.debug(f"Brackets removed. Original: {smi}, New: {new_smi}")
        return new_smi

if __name__ == '__main__':

    cleaner = BracketCleaner(log_level='DEBUG')
    test_smi = '[C]1=[C]C[C@@H](c2ccccc2)N=C1'

    result = cleaner.clean_brackets(test_smi)

    cleaner = BracketCleaner(log_level='DEBUG')
    test_smi = 'Cl[Mg]Cc1ccccc1'

    result = cleaner.clean_brackets(test_smi)