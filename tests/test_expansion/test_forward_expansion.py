from rdkit.Chem import AllChem

from rbc2 import RetroBioCatExpander
from rbc2.configs import Expansion_Config


def test_can_apply_reductive_amination_forwards():

    expansion_config = Expansion_Config()
    expansion_config.use_rdchiral = False  # rdchiral can't handle multiple reactants
    expansion_config.rule_applicator_log_level = 'DEBUG'
    rbc_expander = RetroBioCatExpander(reverse_rules=True,
                                       only_named_reactions=['Reductive amination'],
                                       config=expansion_config)

    smi = 'CCC=O.NC'

    #outcomes = rbc_expander.get_reactions(smi)
    #assert outcomes[0].product == 'CCCNC'

    outcomes = rbc_expander.get_reactions("O=C1CCCCC1.C1CCNCC1")
    print(outcomes)