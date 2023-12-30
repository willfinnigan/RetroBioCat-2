from rbc2.precedent_identification.data_retrieval.bkms.bkms_precedent_data import BKMS_Data
from rbc2.precedent_identification.data_retrieval.enzymemap.enzymemap_precedent_data import EnzymeMap_Data


def test_can_retrieve_data_for_bkms_precedents():
    bkms_data = BKMS_Data()
    result = bkms_data.query_data(query_ids=[1,2,3])
    assert len(result) > 0

def test_can_retrieve_data_for_enzyme_map_precedents():
    em_data = EnzymeMap_Data()
    result = em_data.query_data(template_id=1)
    assert len(result) > 0