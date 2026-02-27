import pytest

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.reaction import Reaction, reaction_from_dict
from rbc2.data_model.reaction_option import ReactionOption, create_evaluate_option_method
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rbc2.template_application.apply_template.rdchiral.main import rdchiralRun
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.components.multi_step_applicator import Multi_step_applicator
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator
from rbc2.template_application.create_reactions_from_output.create_reactions import create_reactions


class TestRdchiralRunAtomMapping:
    """Tests for atom mapping at the rdchiralRun level."""

    def test_single_hit_returns_mapping(self):
        """Simple carbonyl reduction: one product, one mapping."""
        rxn = rdchiralReaction('[C:1]=[O:2]>>[C:1][OH:2]')
        reactants = rdchiralReactants('C=O')  # formaldehyde
        outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

        assert len(outcomes) == 1
        assert len(atom_mappings) == 1

        smi = outcomes[0]
        mapping = atom_mappings[smi]
        # Template has map nums 1 (C) and 2 (O)
        assert 1 in mapping
        assert 2 in mapping
        # Values should be valid atom indices (0-based)
        assert all(isinstance(v, int) and v >= 0 for v in mapping.values())

    def test_diketone_multi_hit_distinct_mappings(self):
        """Diketone O=CC(=O)C has two C=O groups; each outcome should map to a different atom."""
        rxn = rdchiralReaction('[C:1]=[O:2]>>[C:1][OH:2]')
        reactants = rdchiralReactants('O=CC(=O)C')
        outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

        assert len(outcomes) >= 2
        assert len(atom_mappings) >= 2

        # Collect the C atom indices (template mapnum 1) across outcomes
        c_indices = set()
        for smi in outcomes:
            mapping = atom_mappings[smi]
            assert 1 in mapping
            assert 2 in mapping
            c_indices.add(mapping[1])

        # The two outcomes should map to different C atoms
        assert len(c_indices) >= 2

    def test_no_match_returns_empty(self):
        """Template that doesn't match should return empty lists."""
        rxn = rdchiralReaction('[Cl:1]>>[Cl:1]C')
        reactants = rdchiralReactants('C=O')
        outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

        assert outcomes == []
        assert atom_mappings == {}

    def test_backwards_compat_without_flag(self):
        """Without return_atom_mapping, return type should be unchanged (list)."""
        rxn = rdchiralReaction('[C:1]=[O:2]>>[C:1][OH:2]')
        reactants = rdchiralReactants('C=O')
        result = rdchiralRun(rxn, reactants)

        assert isinstance(result, list)

    def test_return_mapped_and_atom_mapping(self):
        """Both return_mapped and return_atom_mapping can be used together."""
        rxn = rdchiralReaction('[C:1]=[O:2]>>[C:1][OH:2]')
        reactants = rdchiralReactants('C=O')
        outcomes, mapped_outcomes, atom_mappings = rdchiralRun(
            rxn, reactants, return_mapped=True, return_atom_mapping=True
        )

        assert isinstance(outcomes, list)
        assert isinstance(mapped_outcomes, dict)
        assert isinstance(atom_mappings, dict)
        assert len(outcomes) == 1
        assert len(atom_mappings) == 1

    def test_mapping_values_are_valid_indices(self):
        """Mapped atom indices should be within range of the target molecule's atoms."""
        rxn = rdchiralReaction('[C:1]=[O:2]>>[C:1][OH:2]')
        reactants = rdchiralReactants('CC=O')  # acetaldehyde, 4 atoms (C, C, =O, implicit H not counted)
        outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

        assert len(outcomes) >= 1
        smi = outcomes[0]
        mapping = atom_mappings[smi]
        # acetaldehyde has 3 heavy atoms (indices 0, 1, 2)
        for idx in mapping.values():
            assert 0 <= idx < 3


# ============================================================
# A. Positional tracking through parser
# ============================================================

class TestParseWithMappings:
    """Tests for _parse_with_mappings in RdChiral_Applicator."""

    def test_parse_with_mappings_basic(self):
        """Simple case — mapping survives split/clean/canonicalize."""
        applicator = RdChiral_Applicator()
        products = ['CO']
        mappings = [{1: 0, 2: 1}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        assert len(parsed) == 1
        assert len(parsed_mappings) == 1
        assert parsed_mappings[0] == {1: 0, 2: 1}

    def test_parse_with_mappings_bracket_cleaning(self):
        """Raw SMILES with unnecessary brackets — mapping still attached after cleaning."""
        config = Expansion_Config()
        config.clean_brackets = True
        applicator = RdChiral_Applicator(config=config)

        products = ['[CH4]']
        mappings = [{1: 0}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        assert len(parsed) == 1
        assert len(parsed_mappings) == 1
        assert parsed_mappings[0] == {1: 0}

    def test_parse_with_mappings_rdkit_canonicalization(self):
        """Non-canonical SMILES, mapping still attached after canonicalization."""
        config = Expansion_Config()
        config.force_rdkit_smis = True
        applicator = RdChiral_Applicator(config=config)

        # OC is non-canonical for CO (methanol)
        products = ['OC']
        mappings = [{1: 0, 2: 1}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        assert len(parsed) == 1
        assert len(parsed_mappings) == 1
        assert parsed_mappings[0] == {1: 0, 2: 1}

    def test_parse_with_mappings_multi_substrate(self):
        """Dot-joined products correctly split, mapping preserved."""
        applicator = RdChiral_Applicator()
        products = ['CO.CC']
        mappings = [{1: 0, 2: 1}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        assert len(parsed) == 1
        assert len(parsed_mappings) == 1
        # products should be split into ['CO', 'CC']
        assert len(parsed[0]) == 2
        assert parsed_mappings[0] == {1: 0, 2: 1}

    def test_parse_with_mappings_dedup(self):
        """Duplicate outcomes — first mapping survives."""
        applicator = RdChiral_Applicator()
        products = ['CO', 'CO']
        mappings = [{1: 0, 2: 1}, {1: 2, 2: 3}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        assert len(parsed) == 1
        assert len(parsed_mappings) == 1
        # First mapping kept
        assert parsed_mappings[0] == {1: 0, 2: 1}

    def test_parse_with_mappings_drop_nones(self):
        """Invalid SMILES dropped, mapping dropped too."""
        config = Expansion_Config()
        config.force_rdkit_smis = True
        applicator = RdChiral_Applicator(config=config)

        # 'invalid_smiles' will become None after rdkit canonicalization
        products = ['CO', 'invalid_smiles']
        mappings = [{1: 0}, {1: 5}]
        parsed, parsed_mappings = applicator._parse_with_mappings(products, mappings)

        # Only the valid one should survive
        assert len(parsed) == len(parsed_mappings)
        assert all(None not in p for p in parsed)


class TestApplicatorAtomMapping:
    """Tests for atom mapping at the RdChiral_Applicator level (parallel lists)."""

    def test_apply_rules_with_atom_mapping(self):
        applicator = RdChiral_Applicator()
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        rxns_dict = {'carbonyl_reduction': rxns}
        smi = 'C=O'

        products, atom_mappings = applicator.apply_rules(smi, rxns_dict, return_atom_mapping=True)

        assert 'carbonyl_reduction' in products
        assert 'carbonyl_reduction' in atom_mappings
        # Parallel: same number of products and mappings
        assert len(products['carbonyl_reduction']) == len(atom_mappings['carbonyl_reduction'])

    def test_apply_rules_without_flag_unchanged(self):
        """Without return_atom_mapping, return is just the product dict."""
        applicator = RdChiral_Applicator()
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        rxns_dict = {'carbonyl_reduction': rxns}
        smi = 'C=O'

        result = applicator.apply_rules(smi, rxns_dict)
        assert isinstance(result, dict)
        assert 'carbonyl_reduction' in result

    def test_apply_reactions_with_atom_mapping(self):
        applicator = RdChiral_Applicator()
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        reactant = applicator.get_reactant('C=O')

        products, atom_mappings = applicator.apply_reactions(reactant, rxns, return_atom_mapping=True)

        assert len(products) > 0
        # Parallel lists
        assert isinstance(atom_mappings, list)
        assert len(atom_mappings) == len(products)

    def test_no_match_returns_empty_with_atom_mapping(self):
        applicator = RdChiral_Applicator()
        rxns = applicator.get_rxns(["[Cl:1]>>[Cl:1]C"])
        rxns_dict = {'chloride': rxns}
        smi = 'C=O'

        products, atom_mappings = applicator.apply_rules(smi, rxns_dict, return_atom_mapping=True)

        assert products == {}
        assert atom_mappings == {}

    def test_diketone_parallel_mappings(self):
        """Diketone with two C=O groups: each outcome gets its own mapping."""
        applicator = RdChiral_Applicator()
        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        rxns_dict = {'carbonyl_reduction': rxns}
        smi = 'O=CC(=O)C'

        products, atom_mappings = applicator.apply_rules(smi, rxns_dict, return_atom_mapping=True)

        assert 'carbonyl_reduction' in products
        product_list = products['carbonyl_reduction']
        mapping_list = atom_mappings['carbonyl_reduction']

        assert len(product_list) == len(mapping_list)
        assert len(product_list) >= 2

        # Each mapping should reference different C atoms
        c_indices = set()
        for m in mapping_list:
            assert 1 in m  # template :1 is C
            c_indices.add(m[1])
        assert len(c_indices) >= 2


class TestRuleApplicatorAtomMapping:
    """Tests for atom mapping at the RuleApplicator level (parallel dicts)."""

    def test_apply_with_atom_mapping(self):
        applicator = RuleApplicator()
        rxns = {'carbonyl_reduction': applicator.smarts_to_rdchiral(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])}
        smi = 'C=O'

        products, atom_mappings = applicator.apply(smi, rxns, return_atom_mapping=True)

        assert 'carbonyl_reduction' in products
        assert 'carbonyl_reduction' in atom_mappings
        assert len(products['carbonyl_reduction']) == len(atom_mappings['carbonyl_reduction'])

    def test_apply_without_flag_unchanged(self):
        """Without return_atom_mapping, return is just the product dict."""
        applicator = RuleApplicator()
        rxns = {'carbonyl_reduction': applicator.smarts_to_rdchiral(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])}
        smi = 'C=O'

        result = applicator.apply(smi, rxns)
        assert isinstance(result, dict)

    def test_template_flag_filtering_removes_mapping(self):
        """Outcome removed by flag check doesn't cause index mismatch."""
        applicator = RuleApplicator()
        # Use a rule that produces a cyclic outcome (should be filtered)
        # and a normal rule
        rxns = {'carbonyl_reduction': applicator.smarts_to_rdchiral(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])}
        smi = 'CC=O'

        products, atom_mappings = applicator.apply(smi, rxns, return_atom_mapping=True)

        # Products and mappings should still be parallel
        for name in products:
            if name in atom_mappings:
                assert len(products[name]) == len(atom_mappings[name])


# ============================================================
# B. create_reactions integration — step-based format
# ============================================================

class TestCreateReactionsAtomMapping:

    def test_create_reactions_stores_atom_mapping_step_format(self):
        """Mapping stored in Reaction.data['atom_mapping'] as list of step dicts."""
        outcomes = {'reduction': [['CO']]}
        atom_mappings = {'reduction': [{1: 0, 2: 1}]}

        reactions = create_reactions('C=O', outcomes, atom_mappings=atom_mappings)

        assert len(reactions) == 1
        steps = reactions[0].data['atom_mapping']
        assert isinstance(steps, list)
        assert len(steps) == 1
        assert steps[0]['rxn_name'] == 'reduction'
        assert steps[0]['target_smi'] == 'C=O'
        assert steps[0]['substrate_smis'] == ['CO']
        assert steps[0]['mapping'] == {1: 0, 2: 1}

    def test_create_reactions_no_mapping_gives_no_key(self):
        """Without atom_mappings, data has no 'atom_mapping' key."""
        outcomes = {'reduction': [['CO']]}

        reactions = create_reactions('C=O', outcomes)

        assert len(reactions) == 1
        assert 'atom_mapping' not in reactions[0].data

    def test_create_reactions_multiple_outcomes_get_correct_mappings(self):
        """Two outcomes with different mappings — each Reaction gets the right one."""
        outcomes = {'reduction': [['CO'], ['OC']]}
        atom_mappings = {'reduction': [{1: 0, 2: 1}, {1: 2, 2: 3}]}

        reactions = create_reactions('C=O', outcomes, atom_mappings=atom_mappings)

        assert len(reactions) == 2
        assert reactions[0].data['atom_mapping'][0]['mapping'] == {1: 0, 2: 1}
        assert reactions[1].data['atom_mapping'][0]['mapping'] == {1: 2, 2: 3}

    def test_create_reactions_multi_step_chain_passthrough(self):
        """Multi-step chain (list of step dicts) passed through unchanged."""
        chain = [
            {
                'rxn_name': 'step1',
                'target_smi': 'A',
                'substrate_smis': ['B'],
                'mapping': {1: 0}
            },
            {
                'rxn_name': 'step2',
                'target_smi': 'B',
                'substrate_smis': ['C', 'D'],
                'mapping': {1: 1, 2: 2}
            }
        ]
        outcomes = {'multi_rxn': [['C', 'D']]}
        atom_mappings = {'multi_rxn': [chain]}

        reactions = create_reactions('A', outcomes, atom_mappings=atom_mappings)

        assert len(reactions) == 1
        steps = reactions[0].data['atom_mapping']
        assert len(steps) == 2
        assert steps[0]['rxn_name'] == 'step1'
        assert steps[1]['rxn_name'] == 'step2'


# ============================================================
# C. Full pipeline (ReactionOption.evaluate) — step-based format
# ============================================================

class TestEvaluateOptionAtomMapping:

    def test_evaluate_option_stores_atom_mapping_step_format(self):
        """Config enabled — Reaction has step-based atom_mapping."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        rule_applicator = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(rule_applicator, config)

        option = ReactionOption(
            target_smi='CC=O',
            name='carbonyl_reduction',
            smarts=["[#6:1]=[O:2]>>[#6:1]-[OH:2]"],
            metadata={},
            rxn_type='test',
            rxn_domain='test',
            score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()
        assert len(reactions) >= 1
        for rxn in reactions:
            steps = rxn.data['atom_mapping']
            assert isinstance(steps, list)
            assert len(steps) == 1
            assert steps[0]['rxn_name'] == 'carbonyl_reduction'
            assert steps[0]['target_smi'] == 'CC=O'
            assert steps[0]['substrate_smis'] == rxn.substrates
            assert len(steps[0]['mapping']) > 0

    def test_evaluate_option_config_disabled_no_mapping(self):
        """Config disabled (default) — Reaction.data has no 'atom_mapping'."""
        config = Expansion_Config()
        assert config.store_atom_mappings is False  # default off
        rule_applicator = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(rule_applicator, config)

        option = ReactionOption(
            target_smi='CC=O',
            name='carbonyl_reduction',
            smarts=["[#6:1]=[O:2]>>[#6:1]-[OH:2]"],
            metadata={},
            rxn_type='test',
            rxn_domain='test',
            score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()
        assert len(reactions) >= 1
        for rxn in reactions:
            assert 'atom_mapping' not in rxn.data

    def test_evaluate_option_diketone_distinct_mappings(self):
        """Two Reactions have different mapping dicts."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        rule_applicator = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(rule_applicator, config)

        option = ReactionOption(
            target_smi='O=CC(=O)C',
            name='carbonyl_reduction',
            smarts=["[#6:1]=[O:2]>>[#6:1]-[OH:2]"],
            metadata={},
            rxn_type='test',
            rxn_domain='test',
            score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()
        assert len(reactions) >= 2

        # At least two distinct mappings
        c_indices = set()
        for r in reactions:
            steps = r.data['atom_mapping']
            assert isinstance(steps, list)
            assert len(steps) == 1
            m = steps[0]['mapping']
            assert 1 in m
            c_indices.add(m[1])
        assert len(c_indices) >= 2


# ============================================================
# D. Edge cases — updated for step-based format
# ============================================================

class TestEnantiomerCombinationWithMapping:
    """Tests that enantiomer combination works with atom mapping enabled."""

    def test_enantiomers_combined_with_mapping_butanone(self):
        """R/S chiral templates on butanone: enantiomers combined to 1 product,
        mapping preserved from first enantiomer."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        config.combine_enantiomers = True  # default, but explicit
        ra = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(ra, config)
        option = ReactionOption(
            target_smi='CCC(=O)C',
            name='chiral_reduction',
            smarts=[
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
            ],
            metadata={},
            rxn_type='test', rxn_domain='test', score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()

        # Enantiomers should be combined into 1 achiral product
        assert len(reactions) == 1
        rxn = reactions[0]
        steps = rxn.data['atom_mapping']
        assert isinstance(steps, list)
        assert len(steps) == 1
        m = steps[0]['mapping']
        assert len(m) > 0
        # All 4 template mapnums should be present
        assert 1 in m and 2 in m and 3 in m and 4 in m

    def test_enantiomers_combined_with_mapping_decalone(self):
        """R/S templates on decalone: enantiomers combined, mapping preserved."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        ra = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(ra, config)
        option = ReactionOption(
            target_smi='O=C1CCCC2CCCCC12',
            name='chiral_reduction',
            smarts=[
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
            ],
            metadata={},
            rxn_type='test', rxn_domain='test', score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()

        assert len(reactions) == 1
        rxn = reactions[0]
        steps = rxn.data['atom_mapping']
        assert isinstance(steps, list)
        assert len(steps) == 1
        # Product should be achiral (no @ in substrates)
        assert '@' not in rxn.substrates[0]

    def test_enantiomers_disabled_gives_both(self):
        """With combine_enantiomers=False, both enantiomers kept with distinct mappings."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        config.combine_enantiomers = False
        ra = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(ra, config)
        option = ReactionOption(
            target_smi='CCC(=O)C',
            name='chiral_reduction',
            smarts=[
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
            ],
            metadata={},
            rxn_type='test', rxn_domain='test', score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()

        # Both enantiomers should be separate reactions
        assert len(reactions) == 2
        m0 = reactions[0].data['atom_mapping'][0]['mapping']
        m1 = reactions[1].data['atom_mapping'][0]['mapping']
        # Mappings should differ (mapnums 3 and 4 are swapped)
        assert m0 != m1

    def test_combine_enantiomers_with_mappings_directly(self):
        """Direct test of _combine_enantiomers_with_mappings on the applicator."""
        applicator = RdChiral_Applicator()

        products = ['CC[C@@H](C)O', 'CC[C@H](C)O']
        mappings = [{1: 2, 2: 3, 3: 4, 4: 1}, {1: 2, 2: 3, 3: 1, 4: 4}]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        # Should combine to 1 achiral product
        assert len(combined_products) == 1
        assert len(combined_mappings) == 1
        assert '@' not in combined_products[0]
        # Mapping should have all 4 mapnums
        assert len(combined_mappings[0]) == 4

    def test_combine_enantiomers_no_pairs_unchanged(self):
        """When there are no enantiomer pairs, products and mappings unchanged."""
        applicator = RdChiral_Applicator()

        products = ['CCO', 'CC=O']
        mappings = [{1: 0}, {1: 1}]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert combined_products == products
        assert combined_mappings == mappings

    def test_combine_enantiomers_mixed_chiral_and_achiral(self):
        """Mix of enantiomer pairs and non-chiral products."""
        applicator = RdChiral_Applicator()

        # Two enantiomers + one achiral product
        products = ['CC[C@@H](C)O', 'CCO', 'CC[C@H](C)O']
        mappings = [{1: 2, 2: 3, 3: 4, 4: 1}, {1: 0, 2: 1}, {1: 2, 2: 3, 3: 1, 4: 4}]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        # Should have 2 products: the achiral CCO + the combined enantiomers
        assert len(combined_products) == 2
        assert len(combined_mappings) == 2


class TestEnantiomerCombinationEdgeCases:
    """Stress tests for _combine_enantiomers_with_mappings matching correctness."""

    def test_diketone_rs_templates_two_positions(self):
        """Diketone with R/S templates: each C=O position produces an R/S pair.
        4 products total → combine to 2 achiral products, each with correct positional mapping."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        config.combine_enantiomers = True
        ra = RuleApplicator(config=config)

        eval_fn = create_evaluate_option_method(ra, config)
        # 2,4-pentanedione: two ketone groups
        option = ReactionOption(
            target_smi='CC(=O)CC(=O)C',
            name='chiral_reduction',
            smarts=[
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
                "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
            ],
            metadata={},
            rxn_type='test', rxn_domain='test', score=1.0,
            evaluation_function=eval_fn,
        )

        reactions = option.evaluate()

        # Should have reactions for each C=O position (enantiomers combined at each)
        assert len(reactions) >= 1
        for rxn in reactions:
            steps = rxn.data.get('atom_mapping', [])
            assert len(steps) == 1
            m = steps[0]['mapping']
            assert len(m) > 0, "Mapping should not be empty after enantiomer combination"
            # All template mapnums should be present
            assert 1 in m and 2 in m

        # If we get 2+ reactions, their mappings should point to different C=O sites
        if len(reactions) >= 2:
            c_sites = set()
            for rxn in reactions:
                c_sites.add(rxn.data['atom_mapping'][0]['mapping'][1])
            assert len(c_sites) >= 2, "Different reactions should map to different C=O sites"

    def test_direct_combine_two_enantiomer_pairs_different_positions(self):
        """Two enantiomer pairs from different positions: each pair maps to different target atoms."""
        applicator = RdChiral_Applicator()

        products = [
            'CC[C@@H](C)O',
            'CCC[C@@H](O)C',
            'CC[C@H](C)O',
            'CCC[C@H](O)C',
        ]
        mappings = [
            {1: 2, 2: 3},
            {1: 5, 2: 6},
            {1: 2, 2: 3},
            {1: 5, 2: 6},
        ]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert len(combined_products) == 2
        assert len(combined_mappings) == 2

        for m in combined_mappings:
            assert len(m) > 0

        sites = {m[1] for m in combined_mappings}
        assert len(sites) == 2, f"Expected 2 distinct sites, got {sites}"

    def test_direct_combine_with_diastereomers(self):
        """Products with two stereocenters."""
        from rbc2.template_application.apply_template.rdchiral.clean import combine_enantiomers_into_racemic

        applicator = RdChiral_Applicator()

        products = [
            'O[C@@H](F)[C@@H](Cl)Br',
            'O[C@H](F)[C@@H](Cl)Br',
            'O[C@@H](F)[C@H](Cl)Br',
            'O[C@H](F)[C@H](Cl)Br',
        ]
        mappings = [
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
        ]

        combined_set = combine_enantiomers_into_racemic(set(products))

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert len(combined_products) == len(combined_mappings)
        for m in combined_mappings:
            assert len(m) > 0, f"Mapping should not be empty, products: {combined_products}"

    def test_direct_combine_potential_ambiguous_stereo_stripped(self):
        """Stress test for stereo-stripped matching."""
        from rbc2.template_application.apply_template.rdchiral.clean import combine_enantiomers_into_racemic

        applicator = RdChiral_Applicator()

        r_smi = '[C@@H](O)(Cl)C(C)C'
        s_smi = '[C@H](O)(Cl)C(C)C'

        combined = combine_enantiomers_into_racemic({r_smi, s_smi})
        if len(combined) == 1:
            products = [r_smi, s_smi]
            mappings = [{1: 0, 2: 1}, {1: 0, 2: 1}]

            combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
                products, mappings
            )

            assert len(combined_products) == 1
            assert len(combined_mappings) == 1
            assert combined_mappings[0] == {1: 0, 2: 1}
        else:
            products = [r_smi, s_smi]
            mappings = [{1: 0, 2: 1}, {1: 0, 2: 1}]

            combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
                products, mappings
            )
            assert len(combined_products) == len(combined_mappings)

    def test_cis_trans_combination_with_mapping(self):
        """E/Z pair combination test."""
        from rbc2.template_application.apply_template.rdchiral.clean import combine_enantiomers_into_racemic

        applicator = RdChiral_Applicator()

        e_smi = 'C/C=C/C'
        z_smi = r'C/C=C\C'

        combined = combine_enantiomers_into_racemic({e_smi, z_smi})

        products = [e_smi, z_smi]
        mappings = [{1: 0, 2: 1}, {1: 0, 2: 1}]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert len(combined_products) == len(combined_mappings)
        for m in combined_mappings:
            assert len(m) > 0

    def test_mixed_tetrahedral_and_ez_products(self):
        """Products where some have tetrahedral stereo and some have E/Z."""
        applicator = RdChiral_Applicator()

        products = [
            'CC[C@@H](C)O',
            'CC[C@H](C)O',
            'C/C=C/C',
        ]
        mappings = [
            {1: 0, 2: 1},
            {1: 0, 2: 1},
            {1: 3, 2: 4},
        ]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert len(combined_products) == 2
        assert len(combined_mappings) == 2

    def test_full_pipeline_asymmetric_diketone_with_rs(self):
        """Asymmetric diketone with R/S templates."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        config.combine_enantiomers = True

        applicator = RdChiral_Applicator(config=config)

        target = 'CC(=O)CCC=O'
        rxns = applicator.get_rxns([
            "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
            "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
        ])
        reactant = applicator.get_reactant(target)

        products, atom_mappings = applicator.apply_reactions(
            reactant, rxns, return_atom_mapping=True
        )

        assert len(products) == len(atom_mappings)
        for m in atom_mappings:
            assert len(m) > 0, "Every product should have a mapping after enantiomer combination"

    def test_asymmetric_diketone_both_positions_mapped_correctly(self):
        """Asymmetric diketone: R/S templates produce enantiomers at BOTH C=O positions."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        config.combine_enantiomers = True

        applicator = RdChiral_Applicator(config=config)
        rxns = applicator.get_rxns([
            '[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]',
            '[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]',
        ])
        reactant = applicator.get_reactant('CCC(=O)CC(=O)C')

        products, atom_mappings = applicator.apply_reactions(
            reactant, rxns, return_atom_mapping=True
        )

        assert len(products) == len(atom_mappings)
        assert len(products) >= 2

        sites = set()
        for m in atom_mappings:
            assert 1 in m and 2 in m
            sites.add(m[1])
        assert len(sites) >= 2, f"Expected 2 distinct C=O sites, got {sites}"

    def test_multi_stereocenter_same_position_safe_mapping(self):
        """When combine_enantiomers partially strips stereo."""
        applicator = RdChiral_Applicator()

        products = [
            'O[C@@H](F)[C@@H](Cl)Br',
            'O[C@H](F)[C@@H](Cl)Br',
            'O[C@@H](F)[C@H](Cl)Br',
            'O[C@H](F)[C@H](Cl)Br',
        ]
        mappings = [
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
            {1: 0, 2: 1, 3: 2, 4: 3},
        ]

        combined_p, combined_m = applicator._combine_enantiomers_with_mappings(products, mappings)

        assert len(combined_p) == 2
        assert len(combined_m) == 2
        for m in combined_m:
            assert m[1] == 0 and m[2] == 1

    def test_combine_enantiomers_empty_mapping_fallback(self):
        """Fallback should be an empty dict (not a crash)."""
        applicator = RdChiral_Applicator()

        products = ['CC[C@@H](C)O', 'CC[C@H](C)O']
        mappings = [{1: 10, 2: 20}, {1: 30, 2: 40}]

        combined_products, combined_mappings = applicator._combine_enantiomers_with_mappings(
            products, mappings
        )

        assert len(combined_products) == 1
        assert len(combined_mappings) == 1
        assert combined_mappings[0] == {1: 10, 2: 20}


class TestAtomMappingEdgeCases:

    def test_atom_mapping_with_stereochemistry(self):
        """Chiral outcomes — products and mappings stay parallel."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        rule_applicator = RuleApplicator(config=config)

        rxns = {'imine_reduction': rule_applicator.smarts_to_rdchiral(
            ["[#6:1]-[#6:2]=[N:3]>>[#6:1]-[#6:2]-[NH:3]"]
        )}

        products, atom_mappings = rule_applicator.apply(
            'c1ccc(C2=NCCCC2)cc1', rxns, return_atom_mapping=True
        )

        if 'imine_reduction' in products:
            product_list = products['imine_reduction']
            mapping_list = atom_mappings['imine_reduction']
            assert len(product_list) == len(mapping_list)

    def test_atom_mapping_serialization_roundtrip(self):
        """Reaction.to_dict() -> reaction_from_dict() preserves atom_mapping in step format."""
        reaction = Reaction(
            product='C=O',
            substrates=['CO'],
            name='test',
            data={'atom_mapping': [
                {
                    'rxn_name': 'test',
                    'target_smi': 'C=O',
                    'substrate_smis': ['CO'],
                    'mapping': {1: 0, 2: 1}
                }
            ]}
        )

        d = reaction.to_dict()
        restored = reaction_from_dict(d)

        steps = restored.data['atom_mapping']
        assert isinstance(steps, list)
        assert len(steps) == 1
        assert steps[0]['rxn_name'] == 'test'
        assert steps[0]['mapping'] == {1: 0, 2: 1}

    def test_symmetric_molecule_one_mapping(self):
        """Symmetric molecule — dedup keeps one Reaction with valid mapping."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        rule_applicator = RuleApplicator(config=config)

        rxns = {'reduction': rule_applicator.smarts_to_rdchiral(
            ["[#6:1](-[#6:2])=[O:3]>>[#6:1](-[#6:2])-[OH:3]"]
        )}

        products, atom_mappings = rule_applicator.apply(
            'CC(=O)C', rxns, return_atom_mapping=True
        )

        if 'reduction' in products:
            assert len(products['reduction']) == len(atom_mappings['reduction'])
            for m in atom_mappings['reduction']:
                assert isinstance(m, dict)
                assert len(m) > 0


# ============================================================
# E. Multi-step atom mapping tests
# ============================================================

class TestMultiStepApplicator:
    """Tests for _apply_steps_with_mappings and apply_multi_step_rules."""

    def _make_applicator_and_rxns(self):
        """Helper: create RdChiral_Applicator and some simple rxns."""
        applicator = RdChiral_Applicator()
        # Imine reduction: C=N -> C-NH
        imine_rxns = applicator.get_rxns([
            '[#6X4;z1:2]-[#7X3;z0:3]>>[#6X3;z1:2]=[#7X2;z0:3]'
        ])
        # Amine oxidation: C-NH -> C=N (not the same — a different step)
        amine_rxns = applicator.get_rxns([
            '[#6:1]-[#6:2]=[N:3]>>[#6:1]-[#6:2]-[NH:3]'
        ])
        return applicator, imine_rxns, amine_rxns

    def test_apply_steps_with_mappings_single_step(self):
        """Single step produces a chain with 1 step dict."""
        applicator = RdChiral_Applicator()
        multi_app = Multi_step_applicator()

        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        step_names = ['carbonyl_reduction']

        products, chains = multi_app._apply_steps_with_mappings(
            applicator, ['CC=O'], [rxns], step_names)

        assert len(products) >= 1
        assert len(chains) == len(products)
        for chain in chains:
            assert len(chain) == 1
            assert chain[0]['rxn_name'] == 'carbonyl_reduction'
            assert chain[0]['target_smi'] == 'CC=O'
            assert isinstance(chain[0]['mapping'], dict)
            assert len(chain[0]['mapping']) > 0

    def test_apply_multi_step_rules_with_mapping_old_format(self):
        """Backward compat: old format (plain list, no tuple) still works."""
        applicator = RdChiral_Applicator()
        multi_app = Multi_step_applicator()

        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])

        # Old format: dict of plain list (no step names)
        multi_step_rxns = {
            'test_rxn': [[rxns]]
        }

        result = multi_app.apply_multi_step_rules(applicator, 'CC=O', multi_step_rxns)
        assert isinstance(result, dict)  # no crash

    def test_apply_multi_step_rules_with_mapping_tuple_format(self):
        """New tuple format returns mappings when return_atom_mapping=True."""
        applicator = RdChiral_Applicator()
        multi_app = Multi_step_applicator()

        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])

        # New format: tuple of (grouped_list_steps, step_names_per_group)
        multi_step_rxns = {
            'test_rxn': ([[rxns]], [['carbonyl_reduction']])
        }

        products, mappings = multi_app.apply_multi_step_rules(
            applicator, 'CC=O', multi_step_rxns, return_atom_mapping=True)

        assert 'test_rxn' in products
        assert 'test_rxn' in mappings
        assert len(products['test_rxn']) == len(mappings['test_rxn'])

        for chain in mappings['test_rxn']:
            assert isinstance(chain, list)
            assert len(chain) == 1
            assert chain[0]['rxn_name'] == 'carbonyl_reduction'

    def test_apply_multi_step_rules_no_mapping_flag(self):
        """Without return_atom_mapping, returns just the product dict (backward compat)."""
        applicator = RdChiral_Applicator()
        multi_app = Multi_step_applicator()

        rxns = applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])

        multi_step_rxns = {
            'test_rxn': ([[rxns]], [['carbonyl_reduction']])
        }

        result = multi_app.apply_multi_step_rules(
            applicator, 'CC=O', multi_step_rxns)

        assert isinstance(result, dict)
        assert 'test_rxn' in result


class TestMultiStepThroughRuleApplicator:
    """Tests for multi-step atom mapping through the full RuleApplicator pipeline."""

    def test_multi_step_via_rule_applicator(self):
        """Multi-step rxns through RuleApplicator.apply() with return_atom_mapping=True."""
        config = Expansion_Config()
        config.store_atom_mappings = True
        ra = RuleApplicator(config=config)

        rxns_single = {'carbonyl_reduction': ra.smarts_to_rdchiral(
            ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"]
        )}

        # Build a multi-step rxn using tuple format (one group, one step)
        step_rxns = ra.smarts_to_rdchiral(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"])
        multistep_rxns = {
            'multi_carbonyl': ([[step_rxns]], [['carbonyl_sub_step']])
        }

        products, atom_mappings = ra.apply(
            'CC=O', rxns_single, multistep_rxns=multistep_rxns, return_atom_mapping=True)

        # Should have both single and multi-step results
        assert 'carbonyl_reduction' in products
        assert 'multi_carbonyl' in products

        # Multi-step mappings should be chains (lists of step dicts)
        if 'multi_carbonyl' in atom_mappings:
            for chain in atom_mappings['multi_carbonyl']:
                assert isinstance(chain, list)
                assert len(chain) >= 1

    def test_multi_step_mixed_with_single_in_create_reactions(self):
        """Single-step gets 1-element chain, multi-step gets N-element chain."""
        # Simulate outcomes from both single and multi-step
        outcomes = {
            'single_rxn': [['CO']],
            'multi_rxn': [['C', 'D']],
        }

        # Single-step mapping is a flat dict
        # Multi-step mapping is a chain (list of step dicts)
        multi_chain = [
            {'rxn_name': 'step1', 'target_smi': 'A', 'substrate_smis': ['B'], 'mapping': {1: 0}},
            {'rxn_name': 'step2', 'target_smi': 'B', 'substrate_smis': ['C', 'D'], 'mapping': {1: 1}},
        ]
        atom_mappings = {
            'single_rxn': [{1: 0, 2: 1}],
            'multi_rxn': [multi_chain],
        }

        reactions = create_reactions('A', outcomes, atom_mappings=atom_mappings)

        # Find the single and multi reactions
        single_rxns = [r for r in reactions if r.name == 'single_rxn']
        multi_rxns = [r for r in reactions if r.name == 'multi_rxn']

        assert len(single_rxns) == 1
        assert len(multi_rxns) == 1

        # Single: 1 step dict with rxn_name from dict key
        single_steps = single_rxns[0].data['atom_mapping']
        assert len(single_steps) == 1
        assert single_steps[0]['rxn_name'] == 'single_rxn'

        # Multi: 2 step dicts passed through
        multi_steps = multi_rxns[0].data['atom_mapping']
        assert len(multi_steps) == 2
        assert multi_steps[0]['rxn_name'] == 'step1'
        assert multi_steps[1]['rxn_name'] == 'step2'
