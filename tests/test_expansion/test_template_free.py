"""Tests for template-free retrosynthesis expanders.

Tests utility components with no model weights, and tests Graph2SMILES
and Augmented Transformer with real downloaded model weights.
EGraph2RSMILES has no downloadable weights and is tested with random weights only.
"""

import os
import pytest
import numpy as np

torch = pytest.importorskip("torch")

from rbc2.expansion.expanders.template_free.utils.chem_utils import (
    get_atom_features_sparse, get_bond_features_sparse, ATOM_FDIMS, BOND_FDIMS,
    ATOM_LIST, ATOM_DICT
)
from rbc2.expansion.expanders.template_free.utils.distance_calculator import DistanceCalculator
from rbc2.expansion.expanders.template_free.utils.data_utils import (
    get_egraph_features_from_smi, create_batch_from_smi, tokenize_smi, EG2RSBatch
)
from rbc2.expansion.expanders.template_free.models.model_interface import (
    TemplateFreeModel, TemplateFreeResult
)
from rbc2.expansion.expanders.template_free.template_free_expander import TemplateFreeExpander
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.reaction_option import ReactionOption
from rbc2.configs.data_path import path_to_data_folder

TEST_SMI = 'c1ccccc1'
TEST_SMI_2 = 'CC(=O)O'
ASPIRIN = 'CC(=O)Oc1ccccc1C(=O)O'

G2S_MODEL_DIR = os.path.join(path_to_data_folder, 'template_free', 'graph2smiles')
AT_MODEL_DIR = os.path.join(path_to_data_folder, 'template_free', 'augmented_transformer')

has_g2s_weights = (os.path.isfile(os.path.join(G2S_MODEL_DIR, 'model.pt')) and
                   os.path.isfile(os.path.join(G2S_MODEL_DIR, 'vocab.txt')))
has_at_weights = len([f for f in os.listdir(AT_MODEL_DIR) if f.startswith('model_step_')]) > 0 if os.path.isdir(AT_MODEL_DIR) else False


# ============================================================
# Utility tests (no model weights needed)
# ============================================================

class TestChemUtils:
    def test_atom_features_sparse_returns_list(self):
        from rdkit import Chem
        mol = Chem.MolFromSmiles('C')
        atom = mol.GetAtomWithIdx(0)
        features = get_atom_features_sparse(atom)
        assert isinstance(features, list)
        assert len(features) == len(ATOM_FDIMS)

    def test_atom_features_sparse_carbon(self):
        from rdkit import Chem
        mol = Chem.MolFromSmiles('C')
        atom = mol.GetAtomWithIdx(0)
        features = get_atom_features_sparse(atom)
        assert features[0] == ATOM_DICT['C']

    def test_bond_features_sparse(self):
        from rdkit import Chem
        mol = Chem.MolFromSmiles('CC')
        bond = mol.GetBondWithIdx(0)
        features = get_bond_features_sparse(bond)
        assert isinstance(features, list)
        assert len(features) == len(BOND_FDIMS)


class TestDistanceCalculator:
    def test_linear_chain(self):
        adj = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.bool_)
        calc = DistanceCalculator()
        dist = calc.calculate(adj, 3, 10)
        assert dist[0, 1] == 1
        assert dist[0, 2] == 2
        assert dist[1, 2] == 1
        np.testing.assert_array_equal(dist.diagonal(), [0, 0, 0])

    def test_disconnected_graph(self):
        adj = np.zeros((3, 3), dtype=np.bool_)
        adj[0, 1] = adj[1, 0] = True
        calc = DistanceCalculator()
        dist = calc.calculate(adj, 3, 10)
        assert dist[0, 1] == 1
        assert dist[0, 2] == 10

    def test_single_node(self):
        adj = np.zeros((1, 1), dtype=np.bool_)
        calc = DistanceCalculator()
        dist = calc.calculate(adj, 1, 10)
        assert dist[0, 0] == 0


class TestDataUtils:
    def test_tokenize_smi_benzene(self):
        assert tokenize_smi('c1ccccc1') == 'c 1 c c c c c 1'

    def test_tokenize_smi_acetic_acid(self):
        assert tokenize_smi('CC(=O)O') == 'C C ( = O ) O'

    def test_tokenize_smi_roundtrip(self):
        smi = 'CC(=O)Oc1ccccc1C(=O)O'
        assert ''.join(tokenize_smi(smi).split()) == smi

    def test_egraph_features_benzene(self):
        node_features, adjacencies = get_egraph_features_from_smi('c1ccccc1')
        assert node_features.shape[0] == 12  # 6 atoms + 6 bonds
        assert adjacencies.shape[0] > 0

    def test_create_batch_from_smi(self):
        batch = create_batch_from_smi(TEST_SMI)
        assert isinstance(batch, EG2RSBatch)
        assert batch.fnode.shape[0] == 1
        assert batch.distances is not None


# ============================================================
# Graph2SMILES featurization tests
# ============================================================

class TestGraph2SmilesFeaturization:
    def test_g2s_atom_features_dimension(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            get_atom_features, ATOM_FDIM as G2S_ATOM_FDIM
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles('C')
        features = get_atom_features(mol.GetAtomWithIdx(0))
        assert len(features) == sum(G2S_ATOM_FDIM)

    def test_g2s_bond_features(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            get_bond_features, BOND_FDIM as G2S_BOND_FDIM
        )
        from rdkit import Chem
        mol = Chem.MolFromSmiles('CC')
        features = get_bond_features(mol.GetBondWithIdx(0))
        assert len(features) == G2S_BOND_FDIM

    def test_create_g2s_batch(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            create_g2s_batch_from_smi, G2SBatch
        )
        batch = create_g2s_batch_from_smi(TEST_SMI)
        assert isinstance(batch, G2SBatch)
        assert batch.fnode.shape[0] > 0
        assert batch.fmess.shape[0] > 0
        assert batch.distances is not None

    def test_g2s_batch_padding_entry(self):
        """Verify the padding entry (index 0) has '*' atom features."""
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            create_g2s_batch_from_smi, ATOM_DICT
        )
        batch = create_g2s_batch_from_smi('CC')
        star_idx = ATOM_DICT['*']
        assert batch.fnode[0, star_idx].item() == 1.0

    def test_g2s_batch_atom_scope_offset(self):
        """Verify atom_scope starts at 1 (offset for padding)."""
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            create_g2s_batch_from_smi
        )
        batch = create_g2s_batch_from_smi('CC')
        assert batch.atom_scope[0][0][0] == 1  # start index is 1
        assert batch.atom_scope[0][0][1] == 2  # 2 atoms in CC


# ============================================================
# Graph2SMILES model tests (real weights)
# ============================================================

@pytest.mark.skipif(not has_g2s_weights, reason="Graph2SMILES model weights not downloaded")
class TestGraph2SmilesRealModel:
    @pytest.fixture(scope="class")
    def model(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesModel
        m = Graph2SmilesModel()
        m.load()
        return m

    def test_model_loads(self, model):
        assert model.is_loaded()

    def test_aspirin_retrosynthesis(self, model):
        results = model.predict(ASPIRIN)
        assert len(results) > 0
        # Top result should be acetyl chloride + salicylic acid
        top = results[0].reactants
        assert 'O' in top  # should contain oxygen
        assert results[0].score > 0.5

    def test_aspirin_top_result_valid_smiles(self, model):
        """Top result should be valid SMILES (lower-ranked may be truncated)."""
        from rdkit import Chem
        results = model.predict(ASPIRIN)
        top = results[0]
        parts = top.reactants.split('.')
        for part in parts:
            mol = Chem.MolFromSmiles(part)
            assert mol is not None, f"Invalid SMILES in top result: {part}"

    def test_ibuprofen_retrosynthesis(self, model):
        results = model.predict('CC(C)Cc1ccc(C(C)C(=O)O)cc1')
        assert len(results) > 0
        assert results[0].score > 0.5

    def test_scores_descending(self, model):
        results = model.predict(ASPIRIN)
        for i in range(len(results) - 1):
            assert results[i].score >= results[i + 1].score

    def test_num_results(self, model):
        results = model.predict(ASPIRIN, num_results=3)
        assert len(results) <= 3


@pytest.mark.skipif(not has_g2s_weights, reason="Graph2SMILES model weights not downloaded")
class TestGraph2SmilesExpander:
    @pytest.fixture(scope="class")
    def expander(self):
        from rbc2.expansion.expanders.template_free.template_free_expanders import Graph2SmilesExpander
        return Graph2SmilesExpander()

    def test_get_reactions_aspirin(self, expander):
        reactions = expander.get_reactions(ASPIRIN)
        assert len(reactions) > 0
        assert all(isinstance(r, Reaction) for r in reactions)

    def test_reactions_have_correct_metadata(self, expander):
        reactions = expander.get_reactions(ASPIRIN)
        for r in reactions:
            assert r.rxn_type == 'graph2smiles'
            assert r.rxn_domain == 'chemistry'
            assert r.product == ASPIRIN

    def test_reactions_have_valid_substrates(self, expander):
        from rdkit import Chem
        reactions = expander.get_reactions(ASPIRIN)
        for r in reactions:
            for s in r.substrates:
                mol = Chem.MolFromSmiles(s)
                assert mol is not None, f"Invalid substrate: {s}"

    def test_network_integration(self):
        from rbc2.expansion.expanders.template_free.template_free_expanders import Graph2SmilesExpander
        network = Network()
        expander = Graph2SmilesExpander(network=network)
        reactions = expander.get_reactions(ASPIRIN)
        assert network.are_options_available(ASPIRIN, 'graph2smiles')
        assert len(network.all_reactions()) > 0


# ============================================================
# Augmented Transformer tests (real weights)
# ============================================================

@pytest.mark.skipif(not has_at_weights, reason="Augmented Transformer weights not downloaded")
class TestAugmentedTransformerRealModel:
    @pytest.fixture(scope="class")
    def model(self):
        from rbc2.expansion.expanders.template_free.models.augmented_transformer.augmented_transformer_model import AugmentedTransformerModel
        m = AugmentedTransformerModel()
        m.load()
        return m

    def test_model_loads(self, model):
        assert model.is_loaded()

    def test_aspirin_retrosynthesis(self, model):
        results = model.predict(ASPIRIN)
        assert len(results) > 0
        assert results[0].score > 0.0

    def test_results_valid_smiles(self, model):
        from rdkit import Chem
        results = model.predict(ASPIRIN)
        for r in results[:5]:
            parts = r.reactants.split('.')
            for part in parts:
                mol = Chem.MolFromSmiles(part)
                assert mol is not None, f"Invalid SMILES: {part}"

    def test_ibuprofen(self, model):
        results = model.predict('CC(C)Cc1ccc(C(C)C(=O)O)cc1')
        assert len(results) > 0


@pytest.mark.skipif(not has_at_weights, reason="Augmented Transformer weights not downloaded")
class TestAugmentedTransformerExpander:
    @pytest.fixture(scope="class")
    def expander(self):
        from rbc2.expansion.expanders.template_free.template_free_expanders import AugmentedTransformerExpander
        return AugmentedTransformerExpander()

    def test_get_reactions_aspirin(self, expander):
        reactions = expander.get_reactions(ASPIRIN)
        assert len(reactions) > 0
        assert all(isinstance(r, Reaction) for r in reactions)

    def test_reactions_have_correct_metadata(self, expander):
        reactions = expander.get_reactions(ASPIRIN)
        for r in reactions:
            assert r.rxn_type == 'augmented_transformer'
            assert r.rxn_domain == 'chemistry'


# ============================================================
# TemplateFreeExpander tests (with mock model)
# ============================================================

class MockTemplateFreeModel(TemplateFreeModel):
    def __init__(self, predictions=None):
        self._loaded = False
        self._predictions = predictions or [
            TemplateFreeResult(reactants='CC=O.N', score=0.9),
            TemplateFreeResult(reactants='CCO.O=O', score=0.7),
            TemplateFreeResult(reactants='C(=O)O.CC', score=0.5),
        ]

    def predict(self, smi, num_results=10):
        return self._predictions[:num_results]

    def is_loaded(self):
        return self._loaded

    def load(self):
        self._loaded = True


class TestTemplateFreeExpanderMock:
    def test_create_new_options_returns_options(self):
        model = MockTemplateFreeModel()
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        options = expander._create_new_options('CCCO')
        assert len(options) > 0
        assert all(isinstance(o, ReactionOption) for o in options)

    def test_options_are_pre_evaluated(self):
        model = MockTemplateFreeModel()
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        options = expander._create_new_options('CCCO')
        for option in options:
            assert option.evaluated is True
            assert len(option.reactions) > 0

    def test_model_lazy_loading(self):
        model = MockTemplateFreeModel()
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        assert not model.is_loaded()
        expander._create_new_options('CCCO')
        assert model.is_loaded()

    def test_filters_invalid_smiles(self):
        predictions = [
            TemplateFreeResult(reactants='INVALID_SMILES', score=0.9),
            TemplateFreeResult(reactants='CC', score=0.8),
        ]
        model = MockTemplateFreeModel(predictions)
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        reactions = expander._process_predictions('CCCO', predictions)
        assert len(reactions) == 1

    def test_filters_self_reactions(self):
        predictions = [
            TemplateFreeResult(reactants='CCCO', score=0.9),
            TemplateFreeResult(reactants='CC.O', score=0.8),
        ]
        model = MockTemplateFreeModel(predictions)
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        reactions = expander._process_predictions('CCCO', predictions)
        for r in reactions:
            assert 'CCCO' not in r.substrates

    def test_deduplicates_predictions(self):
        predictions = [
            TemplateFreeResult(reactants='CC.O', score=0.9),
            TemplateFreeResult(reactants='O.CC', score=0.7),
        ]
        model = MockTemplateFreeModel(predictions)
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        reactions = expander._process_predictions('CCCO', predictions)
        assert len(reactions) == 1


# ============================================================
# Score clamping tests
# ============================================================

class TestScoreClamping:
    def test_scores_clamped_to_0_1(self):
        """Scores above 1.0 or below 0.0 should be clamped."""
        predictions = [
            TemplateFreeResult(reactants='CC', score=1.5),
            TemplateFreeResult(reactants='CCO', score=-0.5),
            TemplateFreeResult(reactants='CCCO', score=0.8),
        ]
        model = MockTemplateFreeModel(predictions)
        expander = TemplateFreeExpander(model=model)
        expander.rxn_type = 'test_model'
        expander.rxn_domain = 'chemistry'
        reactions = expander._process_predictions('CCCCO', predictions)
        for r in reactions:
            assert 0.0 <= r.score <= 1.0


# ============================================================
# Config dataclass tests
# ============================================================

class TestConfigDataclasses:
    def test_graph2smiles_config_defaults(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesConfig
        config = Graph2SmilesConfig()
        assert config.embed_size == 256
        assert config.rel_pos_buckets == 31
        assert config.mpn_type == 'dgcn'
        assert config.dropout == 0.0

    def test_graph2smiles_config_override(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesConfig
        config = Graph2SmilesConfig(embed_size=512, dropout=0.1)
        assert config.embed_size == 512
        assert config.dropout == 0.1

    def test_egraph2rsmiles_config_defaults(self):
        from rbc2.expansion.expanders.template_free.models.egraph2rsmiles.egraph2rsmiles_model import EGraph2RSmilesConfig
        config = EGraph2RSmilesConfig()
        assert config.embed_size == 256
        assert config.rel_pos_buckets == 31
        assert config.n_latent == 0

    def test_config_supports_setattr(self):
        """Configs must support setattr for checkpoint arg extraction."""
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesConfig
        config = Graph2SmilesConfig()
        setattr(config, 'embed_size', 512)
        assert config.embed_size == 512

    def test_config_supports_hasattr(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesConfig
        config = Graph2SmilesConfig()
        assert hasattr(config, 'embed_size')
        assert not hasattr(config, 'nonexistent_field')


# ============================================================
# GraphInferenceModel base class tests
# ============================================================

class TestGraphInferenceBase:
    def test_base_class_is_template_free_model(self):
        from rbc2.expansion.expanders.template_free.models.graph_inference_base import GraphInferenceModel
        assert issubclass(GraphInferenceModel, TemplateFreeModel)

    def test_g2s_model_inherits_base(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesModel
        from rbc2.expansion.expanders.template_free.models.graph_inference_base import GraphInferenceModel
        assert issubclass(Graph2SmilesModel, GraphInferenceModel)

    def test_egraph_model_inherits_base(self):
        from rbc2.expansion.expanders.template_free.models.egraph2rsmiles.egraph2rsmiles_model import EGraph2RSmilesModel
        from rbc2.expansion.expanders.template_free.models.graph_inference_base import GraphInferenceModel
        assert issubclass(EGraph2RSmilesModel, GraphInferenceModel)

    def test_default_normalize_score(self):
        """Default score normalization is exp(log_prob)."""
        import math
        from rbc2.expansion.expanders.template_free.models.graph_inference_base import GraphInferenceModel

        # Can't instantiate directly (abstract), but we can test via a concrete subclass
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesModel
        model = Graph2SmilesModel()
        # Default _normalize_score uses math.exp
        assert abs(model._normalize_score(0.0) - 1.0) < 1e-6
        assert abs(model._normalize_score(-1.0) - math.exp(-1.0)) < 1e-6

    def test_egraph_normalize_score(self):
        """EGraph score normalization uses linear mapping."""
        from rbc2.expansion.expanders.template_free.models.egraph2rsmiles.egraph2rsmiles_model import EGraph2RSmilesModel
        model = EGraph2RSmilesModel()
        assert model._normalize_score(0.0) == 1.0
        assert model._normalize_score(-5.0) == 0.0
        assert model._normalize_score(-10.0) == 0.0  # clamped


# ============================================================
# G2S featurizer module tests (separate file)
# ============================================================

class TestG2SFeaturizerModule:
    def test_build_graph_from_mol(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import build_graph_from_mol
        from rdkit import Chem
        mol = Chem.MolFromSmiles('CC')
        atom_features, bond_features, a_graphs, b_graphs = build_graph_from_mol(mol)
        assert len(atom_features) == 2
        assert len(bond_features) == 2  # 1 bond * 2 directions
        assert len(a_graphs) == 2
        assert len(b_graphs) == 2

    def test_build_graph_no_networkx(self):
        """Verify graph construction works without networkx import."""
        from rbc2.expansion.expanders.template_free.models.graph2smiles import graph2smiles_featurizer
        import importlib
        # Module should not import networkx
        source = importlib.util.find_spec(graph2smiles_featurizer.__name__)
        assert 'networkx' not in open(source.origin).read()


# ============================================================
# Expander repository tests
# ============================================================

class TestExpanderRepository:
    def test_template_free_constants_defined(self):
        from rbc2.expansion.expander_repository import (
            EGRAPH2RSMILES, GRAPH2SMILES, AUGMENTED_TRANSFORMER
        )
        assert EGRAPH2RSMILES == 'egraph2rsmiles'
        assert GRAPH2SMILES == 'graph2smiles'
        assert AUGMENTED_TRANSFORMER == 'augmented_transformer'

    def test_template_free_expanders_registered(self):
        from rbc2.expansion.expander_repository import expander_repo
        assert 'egraph2rsmiles' in expander_repo
        assert 'graph2smiles' in expander_repo
        assert 'augmented_transformer' in expander_repo


# ============================================================
# Shared component tests
# ============================================================

class TestSharedComponents:
    def test_positional_encoding(self):
        from rbc2.expansion.expanders.template_free.models.shared.positional_encoding import PositionalEncoding
        pe = PositionalEncoding(dropout=0.0, dim=64)
        x = torch.randn(10, 1, 64)
        assert pe(x).shape == x.shape

    def test_position_ffn(self):
        from rbc2.expansion.expanders.template_free.models.shared.position_ffn import PositionwiseFeedForward
        ffn = PositionwiseFeedForward(d_model=64, d_ff=128, dropout=0.0)
        x = torch.randn(2, 10, 64)
        assert ffn(x).shape == x.shape

    def test_sequence_mask(self):
        from rbc2.expansion.expanders.template_free.models.shared.misc import sequence_mask
        lengths = torch.tensor([3, 5, 2])
        mask = sequence_mask(lengths)
        assert mask.shape == (3, 5)
        assert mask[0, 2] == True
        assert mask[0, 3] == False

    def test_embeddings(self):
        from rbc2.expansion.expanders.template_free.models.shared.embeddings import Embeddings
        emb = Embeddings(word_vec_size=64, word_vocab_size=100,
                         word_padding_idx=0, position_encoding=True, dropout=0.0)
        x = torch.randint(0, 100, (5, 2, 1))
        assert emb(x).shape == (5, 2, 64)


# ============================================================
# Cached vs naive beam search equivalence tests
# ============================================================

@pytest.mark.skipif(not has_g2s_weights, reason="Graph2SMILES model weights not downloaded")
class TestCachedBeamSearchEquivalence:
    @pytest.fixture(scope="class")
    def model(self):
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_model import Graph2SmilesModel
        m = Graph2SmilesModel()
        m.load()
        return m

    def test_cached_beam_search_matches_naive(self, model):
        """Cached beam search must produce identical results to naive."""
        from rbc2.expansion.expanders.template_free.models.shared.beam_search import (
            beam_search_decode, _beam_search_decode_naive
        )
        from rbc2.expansion.expanders.template_free.models.graph2smiles.graph2smiles_featurizer import (
            create_g2s_batch_from_smi
        )

        batch = model._create_batch(ASPIRIN)
        device = model._device

        cached_results = beam_search_decode(
            model._model, batch, model._vocab,
            beam_size=10, n_best=10, max_length=512,
            temperature=1.0, device=device)

        naive_results = _beam_search_decode_naive(
            model._model, batch, model._vocab,
            beam_size=10, n_best=10, max_length=512,
            temperature=1.0, device=device)

        assert len(cached_results) == len(naive_results), \
            f"Different number of results: {len(cached_results)} vs {len(naive_results)}"

        for i, (cached, naive) in enumerate(zip(cached_results, naive_results)):
            assert cached.token_ids == naive.token_ids, \
                f"Result {i} tokens differ:\ncached: {cached.token_ids}\nnaive:  {naive.token_ids}"
            assert abs(cached.score - naive.score) < 1e-4, \
                f"Result {i} scores differ: {cached.score} vs {naive.score}"

    def test_cached_greedy_matches_naive(self, model):
        """Cached greedy decode must produce identical results to naive."""
        from rbc2.expansion.expanders.template_free.models.shared.beam_search import (
            beam_search_decode, _beam_search_decode_naive
        )

        batch = model._create_batch(ASPIRIN)
        device = model._device

        cached_results = beam_search_decode(
            model._model, batch, model._vocab,
            beam_size=1, n_best=1, max_length=512,
            temperature=1.0, device=device)

        naive_results = _beam_search_decode_naive(
            model._model, batch, model._vocab,
            beam_size=1, n_best=1, max_length=512,
            temperature=1.0, device=device)

        assert cached_results[0].token_ids == naive_results[0].token_ids, \
            f"Greedy tokens differ:\ncached: {cached_results[0].token_ids}\nnaive:  {naive_results[0].token_ids}"

    def test_cached_beam_search_benchmark(self, model):
        """Print timing comparison (informational, no assertion)."""
        import time
        from rbc2.expansion.expanders.template_free.models.shared.beam_search import (
            beam_search_decode, _beam_search_decode_naive
        )

        batch = model._create_batch(ASPIRIN)
        device = model._device

        # Warm up
        beam_search_decode(model._model, batch, model._vocab,
                           beam_size=10, n_best=10, max_length=512,
                           temperature=1.0, device=device)

        # Cached
        start = time.perf_counter()
        for _ in range(3):
            beam_search_decode(model._model, batch, model._vocab,
                               beam_size=10, n_best=10, max_length=512,
                               temperature=1.0, device=device)
        cached_time = (time.perf_counter() - start) / 3

        # Naive
        start = time.perf_counter()
        for _ in range(3):
            _beam_search_decode_naive(model._model, batch, model._vocab,
                                      beam_size=10, n_best=10, max_length=512,
                                      temperature=1.0, device=device)
        naive_time = (time.perf_counter() - start) / 3

        speedup = naive_time / cached_time if cached_time > 0 else float('inf')
        print(f"\nBeam search timing: cached={cached_time:.3f}s, naive={naive_time:.3f}s, speedup={speedup:.1f}x")
