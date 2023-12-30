from typing import Optional, List

import pandas as pd
from rdkit import Chem

from rbc2.configs.data_path import path_to_data_folder
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expanders.action_getters.bkms.bkms_action_getter import BKMS_Action_Getter
from rbc2.expansion.expanders.action_getters.enzymemap.enzymemap_action_getter import \
    EnzymeMap_Action_Getter
from rbc2.expansion.expanders.action_getters.retrorules.retrorules_getter import RetroRules_Getter
from rbc2.expansion.default_expander_interface import DefaultExpander
from rbc2.precedent_identification.data_retrieval.bkms.bkms_precedent_data import BKMS_Data
from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentData
from rbc2.precedent_identification.data_retrieval.enzymemap.enzymemap_precedent_data import EnzymeMap_Data
from rbc2.precedent_identification.similarity_scorer import SimilarityScorer
from rbc2.reaction_network_entities.network import Network
from rbc2.reaction_network_entities.reaction import Reaction


class RetroRulesExpander(DefaultExpander):
    cofactors_df = None
    cofactor_inchis = None

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 rank_by='combined_score',  # similarity, score, combined_score
                 diameter=6,
                 similarity_threshold=0.2,
                 score_threshold=0.2,
                 combined_score_threshold=0.2,
                 max_reactions=100  # max number of similar substrates to consider reactions for
                 ):

        super().__init__(network=network, config=config)

        self.action_getter = RetroRules_Getter(diameter=diameter,
                                               similarity_threshold=similarity_threshold,
                                               score_threshold=score_threshold,
                                               combined_score_threshold=combined_score_threshold,
                                               max_reactions=max_reactions,
                                               log_level='WARNING')
        self.rxn_type = 'retrorules'
        self.rxn_domain = 'biosynthesis'
        self.score_key = rank_by

    @classmethod
    def _load_cofactor_df(cls):
        if cls.cofactors_df is None:
            cls.cofactors_df = pd.read_csv(f'{path_to_data_folder}/retrorules/retrorules_cofactors.tsv',
                                           sep='\t',
                                           names=['inchi', 'name', 'reactions'])
            cls.cofactor_inchis = set(cls.cofactors_df['inchi'])

    def reaction_processing_function(self, reactions: List[Reaction]) -> List[Reaction]:
        """
        This function will be applied when a ReactionOption is evaluated.
        Remove cofactors from reaction substrates.  Save these as cofactors in template metadata """

        self._load_cofactor_df()  # first time this is called will load the df

        for reaction in reactions:
            # get substrates as inchis (as this is how the retorules cofactors are stored)
            substrate_inchis = [Chem.MolToInchi(Chem.MolFromSmiles(smi)) for smi in reaction.substrates]
            inchi_smi_dict = {inchi: smi for inchi, smi in zip(substrate_inchis, reaction.substrates)}

            # get any substrates which are marked as cofactors
            matches = self.cofactor_inchis.intersection(substrate_inchis)

            # remove matches from reaction.substrates
            reaction.substrates = [inchi_smi_dict[inchi] for inchi in inchi_smi_dict if inchi not in matches]

            # get inchi name from self.cofactors_df and save to template metadata
            for inchi in matches:
                name = self.cofactors_df[self.cofactors_df['inchi'] == inchi]['name'].values[0]
                if 'cofactors' not in reaction.template_metadata[reaction.name]:
                    reaction.template_metadata[reaction.name]['cofactors'] = []
                reaction.template_metadata[reaction.name]['cofactors'].append(name)

        return reactions


class BKMSExpander(DefaultExpander):

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 precedent_data: Optional[PrecedentData] = None,
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 allow_multi_product_templates=False,
                 enable_precedent_search=True,
                 similarity_cutoff=0.1
                 ):
        super().__init__(network=network, config=config)
        self.similarity_cutoff = similarity_cutoff
        self.enable_precedent_search = enable_precedent_search
        self.action_getter = BKMS_Action_Getter(cutoff_number=cutoff_number,
                                                cutoff_cumulative=cutoff_cumulative,
                                                allow_multi_product_templates=allow_multi_product_templates)
        self.rxn_type = 'bkms'
        self.rxn_domain = 'biosynthesis'
        self.score_key = 'score'

        if precedent_data is None:
            precedent_data = BKMS_Data()
        self.precedent_scorer = SimilarityScorer(precedent_data)

    def precedent_evaluation_function(self, reaction: Reaction):
        if self.enable_precedent_search == False:
            return None

        # get the reference ids from the templates - to query the precedent dataset
        query_ids = []
        for template_id, template_metadata in reaction.template_metadata.items():
            for ref_id in template_metadata['references']:
                query_ids.append(ref_id)

        reaction.precedents = self.precedent_scorer.get_precedents(target_smi=reaction.product,
                                                                   topn=1,
                                                                   cutoff=self.similarity_cutoff,
                                                                   query_ids=query_ids)


class EnzymeMapExpander(DefaultExpander):
    cofactor_names = None
    cofactor_smis = None

    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 precedent_data: Optional[PrecedentData] = None,
                 cutoff_cumulative=0.995,
                 cutoff_number=50,
                 enable_precedent_search=True,
                 similarity_cutoff=0.1):

        super().__init__(network=network, config=config)
        self.similarity_cutoff = similarity_cutoff
        self.enable_precedent_search = enable_precedent_search
        self.action_getter = EnzymeMap_Action_Getter(cutoff_cumulative=cutoff_cumulative,
                                                     cutoff_number=cutoff_number)
        self.rxn_type = 'enzymemap'
        self.rxn_domain = 'biosynthesis'
        self.score_key = 'score'
        if precedent_data is None:
            precedent_data = EnzymeMap_Data()
        self.precedent_scorer = SimilarityScorer(precedent_data)

        self._load_cofactor_csv()

    @classmethod
    def _load_cofactor_csv(cls):
        if cls.cofactor_names is None:
            df = pd.read_csv(f'{path_to_data_folder}/enzymemap/enzyme_map_cofactors.csv', header=0)
            cls.cofactor_names = pd.Series(df.name.values, index=df.smiles).to_dict()  # cofactor_csv to a dictionary
            cls.cofactor_smis = set(cls.cofactor_names.keys())

    def precedent_evaluation_function(self, reaction: Reaction):
        if self.enable_precedent_search == False:
            return None

        # use the template_id to query the precedent dataset
        ids_to_query = []
        for template_id, template_metadata in reaction.template_metadata.items():
            t_id = template_metadata['template_id']
            ids_to_query.append(t_id)

        ids_to_query = list(set(ids_to_query))

        precedents = []
        for t_id in ids_to_query:
            precedents += self.precedent_scorer.get_precedents(target_smi=reaction.product,
                                                               topn=1,
                                                               cutoff=self.similarity_cutoff,
                                                               template_id=t_id)
        reaction.precedents = precedents

    def reaction_processing_function(self, reactions: List[Reaction]) -> List[Reaction]:
        """
        This function will be applied when a ReactionOption is evaluated.
        Remove cofactors from reaction substrates.  Save these as cofactors in template metadata """

        for reaction in reactions:
            # get any substrates which are marked as cofactors
            matches = self.cofactor_smis.intersection(reaction.substrates)

            # remove matches from reaction.substrates
            reaction.substrates = [smi for smi in reaction.substrates if smi not in matches]

            # save matches to template metadata
            reaction.template_metadata[reaction.name]['cofactors'] = [self.cofactor_names[smi] for smi in matches]

        return reactions


if __name__ == '__main__':
    expander = RetroRulesExpander()
    expander.config.enable_precedent_search = False
    options = expander.get_options('CC(O)[C@@H](O)c1ccccc1')
    print(len(options))
    reactions = expander.get_reactions('CC(O)[C@@H](O)c1ccccc1')
    print(reactions)
    print(len(reactions))

    # expander = BKMSExpander()
    # reactions = expander.get_reactions('CCCO')
    # print(reactions)

    # expander = EnzymeMapExpander()
    # reactions = expander.get_reactions('CCCCCCC=O')
    # print(reactions[0].precedents)

    # for rxn in reactions:
    # print(rxn)
