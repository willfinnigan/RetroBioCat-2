import pytest

from rbc2.precedent_identification.data_retrieval.bkms.bkms_precedent_data import BKMS_DataQuery
from rbc2.precedent_identification.data_retrieval.enzymemap.enzymemap_precedent_data import EnzymeMap_DataQuery
from rbc2.precedent_identification.data_retrieval.retrobiocat.local_data_query import RetroBioCatLocalPrecedentDataQuery
from rbc2.precedent_identification.similarity_scorer import PandasSimilarityScorer
from rbc2.data_model.precedents import Precedent


def test_enzymemap_similarity_scorer_returns_result_with_similarity_to_target_product():
    scorer = PandasSimilarityScorer(EnzymeMap_DataQuery())
    result = scorer.score_data(target_smi='CCCC=O', topn=1, cutoff=0.0, template_id=1)
    assert 'similarity' in result.iloc[0]

def test_bkms_similarity_scorer_returns_result_with_similarity_to_target_product():
    scorer = PandasSimilarityScorer(BKMS_DataQuery())
    result = scorer.score_data(target_smi='CCCC=O', topn=1, cutoff=0.0, query_ids=[1,2,3])
    assert 'similarity' in result.iloc[0]

def test_bkms_similarity_scorer_without_specifying_ids_raises_an_exception():
    scorer = PandasSimilarityScorer(BKMS_DataQuery())
    with pytest.raises(Exception):
        result = scorer.score_data(target_smi='CCCC=O', topn=1, cutoff=0.0)

def test_retrobiocat_similarity_scorer_returns_result_with_similarity_to_target_product():
    scorer = PandasSimilarityScorer(RetroBioCatLocalPrecedentDataQuery())
    result = scorer.score_data(target_smi='CCCC=O', topn=1, cutoff=0.0, reaction_name='Carboxylic acid reduction', enzyme_types=['CAR'])
    assert 'similarity' in result.iloc[0]


def test_can_create_precedents():
    scorer = PandasSimilarityScorer(RetroBioCatLocalPrecedentDataQuery())

    precedents = scorer.get_precedents(target_smi='CCCC=O',
                                       topn=1,
                                       cutoff=0.0,
                                       reaction_name='Carboxylic acid reduction',
                                       enzyme_types=['CAR'])

    assert isinstance(precedents[0], Precedent) == True
