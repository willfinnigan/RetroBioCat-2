import pytest

from rbc2.reaction_evaluation.feasability import default_filter_repo, AIZYNTHFINDER_FILTER
from rbc2.data_model.reaction import Reaction


@pytest.mark.parametrize(["product", "substrates", "expected"],
                         [["COC(=O)c1ccc(C2OC(CO)C(O)C(OCC(=O)O)C2O)cc1", ["COC(=O)c1ccc(C2OC(CO)C(O)C(O)C2O)cc1", 'O=C(O)CCl'], True],
                          ['COC(=O)c1ccc(cc1)C1OC(COCC2=CC=CC=C2)C(O)C(OCC(O)=O)C1O', ['O=C(O)COC1C(O)C(COCc2ccccc2)OC(c2ccc(C(=O)O)cc2)C1O', 'C=[N+]=[N-]'], True],
                          ["O=C(CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F", ["FC(F)(F)c1nnc2n1CCNC2", "O=C(Cl)CC(=O)Cc1cc(F)c(F)cc1F"], True],
                          ["O=CC1=CC=CC=C1", ["NCCCl"], False]])
def test_aizynthfinder_feasability_score(product, substrates, expected):
    reaction = Reaction(product, substrates)
    feasability = default_filter_repo[AIZYNTHFINDER_FILTER](reaction)
    result = feasability > 0.1
    assert result == expected