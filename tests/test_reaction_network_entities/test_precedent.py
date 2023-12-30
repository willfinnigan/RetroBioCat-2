from dataclasses import asdict

from rbc2.reaction_network_entities.precedents import Precedent


def test_save_and_load_precedent():
    precedent = Precedent(name='test',
                          precedent_id='test_id',
                          data={'test': 'test'},
                          similarity=0.5,
                          rxn_smi='CCCO>>CCC=O')

    precedent_dict = asdict(precedent)
    loaded_precedent = Precedent(**precedent_dict)
    assert loaded_precedent.name == precedent.name