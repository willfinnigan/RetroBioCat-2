from typing import List

from rdkit import Chem

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.download_data_files.download_aizynthfinder_filter import does_aizynthfinder_filter_exist, \
    download_aizynthfinder_filter
from rbc2.utils import fingerprints, load_keras_models

data_folder = f'{path_to_data_folder}/aizynthfinder'

class AIZynthFilterModel():

    def __init__(self):
        self.filter_model = None

        if does_aizynthfinder_filter_exist() == False:
            download_aizynthfinder_filter()

    def load_model(self):
        if self.filter_model == None:
            filter_path = data_folder + '/filter_policy_all.hdf5'
            self.filter_model = load_keras_models.LocalKerasModel(filter_path)

    def run_model(self, product_smi: str, substrate_smis: List[str]) -> float:
        if self.filter_model == None:
            self.load_model()

        product_mol = Chem.MolFromSmiles(product_smi)
        substrate_mols = [Chem.MolFromSmiles(smi) for smi in substrate_smis]

        prod_fp = fingerprints.get_mol_fingerprint(product_mol, radius=2, nBits=len(self.filter_model))
        prod_fp = prod_fp.reshape([1, 2048])

        reaction_fp = fingerprints.get_reaction_fingerprint(product_mol, substrate_mols, radius=2,
                                                            nBits=len(self.filter_model))
        reaction_fp = reaction_fp.reshape([1, 2048])

        kwargs = {"input_1": prod_fp, "input_2": reaction_fp}
        feasability = float(self.filter_model.predict(prod_fp, reaction_fp, **kwargs)[0][0])

        return feasability





if __name__ == '__main__':

    product = "COC(=O)c1ccc(C2OC(CO)C(O)C(OCC(=O)O)C2O)cc1"
    substrates = ["COC(=O)c1ccc(C2OC(CO)C(O)C(O)C2O)cc1", 'O=C(O)CCl']
    f = AIZynthFilterModel().run_model(product, substrates)
    print(f"Should be infeasible - {f}")

    product = "COC(=O)c1ccc(C2OC(CO)C(O)C(OCC(=O)O)C2O)cc1"
    substrates = ['COC(=O)c1ccc(cc1)C1OC(COCC2=CC=CC=C2)C(O)C(OCC(O)=O)C1O']
    f = AIZynthFilterModel().run_model(product, substrates)
    print(f"Should be feasible - {f}")


    product = 'COC(=O)c1ccc(cc1)C1OC(COCC2=CC=CC=C2)C(O)C(OCC(O)=O)C1O'
    substrates = ['O=C(O)COC1C(O)C(COCc2ccccc2)OC(c2ccc(C(=O)O)cc2)C1O', 'C=[N+]=[N-]']
    f = AIZynthFilterModel().run_model(product, substrates)
    print(f"Should be feasible - {f}")


    product = "O=C(CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F"
    substrates = ["FC(F)(F)c1nnc2n1CCNC2", "O=C(Cl)CC(=O)Cc1cc(F)c(F)cc1F"]
    f = AIZynthFilterModel().run_model(product, substrates)
    print(f"Should be feasible - {f}")

    product = "O=CC1=CC=CC=C1"
    substrates = ["NCCCl"]
    f = AIZynthFilterModel().run_model(product, substrates)
    print(f"Should be infeasible - {f}")
