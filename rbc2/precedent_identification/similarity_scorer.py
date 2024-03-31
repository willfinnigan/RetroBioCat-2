from typing import List, Optional, Callable, Tuple
import pandas as pd
from rdkit.DataStructs import ExplicitBitVect

from rbc2.precedent_identification.data_retrieval.data_interface import PrecedentData
from rbc2.precedent_identification.similarity_tools import get_single_fp, bulk_similarity
from rbc2.data_model.precedents import make_precedents

Smi = str
DataQueryFunction = Callable[..., pd.DataFrame]
FpQueryFunction = Callable[[List[Smi]], Tuple[List[ExplicitBitVect], List[Smi]]]
RankingFunction = Callable[[pd.DataFrame, int], pd.DataFrame]

class SimilarityScorer():

    def __init__(self,
                 precedent_data: PrecedentData,
                 ranking_function: Optional[RankingFunction] = None):

        self.precedent_data = precedent_data
        self.ranking_function = ranking_function  # optionally ranks results

    def score_data(self, target_smi: Smi, topn: int, cutoff: float, **key_atts):
        """
        Finds similar data to target smi and returns it as a dataframe
        :param target_smi: the smiles to compare with
        :param topn: the number of top enzymes to return per similar molecule
        :param cutoff: the similarity cutoff
        """
        data_df = self.precedent_data.query_data(**key_atts)  # get the relevant data
        target_fp = get_single_fp(target_smi)  # generate fp for target_smi
        smis = list(data_df[self.precedent_data.product_column].unique())  # the smis to compare target_smi with
        smi_sims = bulk_similarity(target_fp, smis, self.precedent_data.get_fps)  # get similarities to smis
        smi_sims = {smi: sim for smi, sim in smi_sims.items() if sim >= cutoff}  # get dict of smi: sim above cutoff

        result_df = data_df[data_df[self.precedent_data.product_column].isin(smi_sims.keys())].copy()  # get df of similiar smis
        result_df['similarity'] = result_df[self.precedent_data.product_column].map(smi_sims)  # add similarity scores
        result_df['similarity'] = result_df['similarity'].apply(lambda x: round(x, 3))  # round similarity to 3 decimal places
        result_df.sort_values('similarity', ascending=False, inplace=True)  # put most similar first
        if self.ranking_function is not None:
            result_df = self.ranking_function(result_df, topn)  # rank and take only top_n rows per smiles

        return result_df

    def get_precedents(self, target_smi: str, cutoff: float, topn: int = 1, **key_atts):
        result_df = self.score_data(target_smi, topn, cutoff, **key_atts)
        precedents = make_precedents(result_df,
                                     self.precedent_data.enzyme_column,
                                     self.precedent_data.id_column)
        return sorted(precedents, key=lambda x: x.similarity, reverse=True)
