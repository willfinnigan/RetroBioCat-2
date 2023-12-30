from rbc2.template_application.apply_template.components.result_parser import Result_Parser


class Test_result_parser():
    parser = Result_Parser(log_level='DEBUG')

    def test_identical_results_are_combined(self):
        products = ['CO', 'CO']
        parsed_products = self.parser.parse(products)
        assert parsed_products == [['CO']]

    def test_joined_molecules_are_split(self):
        products = ['CO.C']
        parsed_products = self.parser.parse(products)
        assert 'CO' in parsed_products[0]
        assert 'C' in parsed_products[0]

    def test_bracket_cleaning(self):
        products = ['[C]1=[C]C[C@@H](c2ccccc2)N=C1', 'C1=CC[C@@H](c2ccccc2)N=C1']
        parsed_products = self.parser.parse(products)
        assert parsed_products == [['C1=CC[C@@H](c2ccccc2)N=C1']]

    def test_bracket_cleaning_with_mg_compound(self):
        products = ["Cl[Mg]Cc1ccccc1"]
        parsed_products = self.parser.parse(products)
        assert parsed_products == [["Cl[Mg]Cc1ccccc1"]]

    def test_that_any_non_rdkit_accepted_smis_are_dropped(self):
        products = ["not_a_smiles", "CCO"]
        parsed_products = self.parser.parse(products)
        assert parsed_products == [['CCO']]




