from typing import Optional

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.network import Network
from rbc2.expansion.expanders.template_free.template_free_expander import TemplateFreeExpander
from rbc2.expansion.expanders.template_free.models.egraph2rsmiles.egraph2rsmiles_model import EGraph2RSmilesModel
from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesModel
from rbc2.expansion.expanders.template_free.models.augmented_transformer.augmented_transformer_model import AugmentedTransformerModel


class EGraph2RSmilesExpander(TemplateFreeExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 num_results: int = 10):
        model = EGraph2RSmilesModel()
        super().__init__(model=model, network=network, config=config,
                         num_results=num_results)
        self.rxn_type = 'egraph2rsmiles'
        self.rxn_domain = 'chemistry'


class Graph2SmilesExpander(TemplateFreeExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 num_results: int = 10):
        model = Graph2SmilesModel()
        super().__init__(model=model, network=network, config=config,
                         num_results=num_results)
        self.rxn_type = 'graph2smiles'
        self.rxn_domain = 'chemistry'


class AugmentedTransformerExpander(TemplateFreeExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 num_results: int = 10):
        model = AugmentedTransformerModel()
        super().__init__(model=model, network=network, config=config,
                         num_results=num_results)
        self.rxn_type = 'augmented_transformer'
        self.rxn_domain = 'chemistry'
