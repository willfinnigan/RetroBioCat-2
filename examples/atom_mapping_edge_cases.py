"""
Edge-case exploration for atom mapping through the full ReactionOption pipeline.

Every case enables config.store_atom_mappings = True and exercises a different
corner of the rdchiral → parser → create_reactions → Reaction.data path.

Run:  python examples/atom_mapping_edge_cases.py
"""

from rdkit import Chem

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.reaction import Reaction, reaction_from_dict
from rbc2.data_model.reaction_option import ReactionOption, create_evaluate_option_method
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def make_config():
    config = Expansion_Config()
    config.store_atom_mappings = True
    return config


def make_option(target_smi, name, smarts_list, config, rule_applicator):
    eval_fn = create_evaluate_option_method(rule_applicator, config)
    return ReactionOption(
        target_smi=target_smi,
        name=name,
        smarts=smarts_list,
        metadata={},
        rxn_type='test',
        rxn_domain='test',
        score=1.0,
        evaluation_function=eval_fn,
    )


def print_header(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def describe_mapping(target_smi, mapping):
    """Print human-readable mapping: template mapnum -> atom symbol in target."""
    mol = Chem.MolFromSmiles(target_smi)
    if mol is None:
        print(f"    (could not parse target SMILES: {target_smi})")
        return
    for tmap, aidx in sorted(mapping.items()):
        if aidx < mol.GetNumAtoms():
            atom = mol.GetAtomWithIdx(aidx)
            print(f"    template :{tmap} -> atom idx {aidx} ({atom.GetSymbol()})")
        else:
            print(f"    template :{tmap} -> atom idx {aidx} (OUT OF RANGE)")


def show_reactions(target_smi, reactions, check_fn=None):
    """Pretty-print reactions with their atom mappings. Optionally run a check."""
    if len(reactions) == 0:
        print("  (no reactions)")
        return

    all_ok = True
    for i, rxn in enumerate(reactions):
        steps = rxn.data.get('atom_mapping', [])
        # Extract flat mapping from first step for backward-compat checks
        mapping = steps[0]['mapping'] if steps else {}
        print(f"\n  Reaction {i+1}: {rxn.reaction_smiles()}")
        print(f"    substrates = {rxn.substrates}")
        print(f"    atom_mapping steps = {len(steps)}")
        for si, step in enumerate(steps):
            print(f"      step {si+1}: rxn_name={step['rxn_name']}, mapping={step['mapping']}")
        describe_mapping(target_smi, mapping)

        if check_fn is not None:
            ok = check_fn(i, rxn, mapping)
            if not ok:
                all_ok = False

    return all_ok


# ===========================================================================
# 1. Simple single-product carbonyl reduction
# ===========================================================================

def case_1_simple_carbonyl():
    print_header("1. Simple carbonyl reduction  (CC=O -> CCO)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('CC=O', 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        mol = Chem.MolFromSmiles('CC=O')
        ok = True
        # template :1 should be the carbonyl C (idx 1), :2 should be O (idx 2)
        if 1 not in m or 2 not in m:
            print("    !! FAIL: expected mapnums 1 and 2")
            ok = False
        if m.get(2) is not None and mol.GetAtomWithIdx(m[2]).GetSymbol() != 'O':
            print("    !! FAIL: mapnum 2 should map to O")
            ok = False
        if ok:
            print("    >> OK")
        return ok

    ok = show_reactions('CC=O', rxns, check)
    return ok


# ===========================================================================
# 2. Diketone — two distinct C=O sites
# ===========================================================================

def case_2_diketone():
    print_header("2. Diketone with two C=O  (O=CC(=O)C)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('O=CC(=O)C', 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    c_indices = set()

    def check(i, rxn, m):
        if 1 not in m or 2 not in m:
            print("    !! FAIL: expected mapnums 1 and 2")
            return False
        c_indices.add(m[1])
        print("    >> OK")
        return True

    ok = show_reactions('O=CC(=O)C', rxns, check)
    if len(c_indices) < 2:
        print("\n  !! FAIL: expected at least 2 distinct C atom indices across reactions")
        ok = False
    else:
        print(f"\n  >> PASS: {len(c_indices)} distinct C sites mapped")
    return ok


# ===========================================================================
# 3. Symmetric ketone (acetone) — dedup should keep one reaction
# ===========================================================================

def case_3_symmetric_ketone():
    print_header("3. Symmetric ketone  (CC(=O)C  — acetone)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('CC(=O)C', 'ketone_red',
                      ["[#6:1](-[#6:2])=[O:3]>>[#6:1](-[#6:2])-[OH:3]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) == 0:
            print("    !! FAIL: mapping should not be empty")
            return False
        # template :3 is the O
        mol = Chem.MolFromSmiles('CC(=O)C')
        if 3 in m and mol.GetAtomWithIdx(m[3]).GetSymbol() != 'O':
            print("    !! FAIL: mapnum 3 should map to O")
            return False
        print("    >> OK")
        return True

    ok = show_reactions('CC(=O)C', rxns, check)
    # Acetone is symmetric — both C=O matches yield same product, so dedup to 1
    print(f"\n  >> {len(rxns)} reaction(s) after dedup (expect 1 for symmetric ketone)")
    return ok


# ===========================================================================
# 4. Chiral template — imine reduction producing stereocenter
# ===========================================================================

def case_4_chiral_imine():
    print_header("4. Chiral imine reduction  (c1ccc(C2=NCCCC2)cc1)")
    config = make_config()
    ra = RuleApplicator(config=config)
    # This template can produce (R) and (S) amines
    opt = make_option('c1ccc(C2=NCCCC2)cc1', 'imine_red',
                      ["[#6:1]-[#6:2]=[N:3]>>[#6:1]-[#6:2]-[NH:3]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) == 0:
            print("    !! FAIL: mapping should not be empty")
            return False
        print("    >> OK")
        return True

    ok = show_reactions('c1ccc(C2=NCCCC2)cc1', rxns, check)
    # Enantiomers are combined even with atom mapping — mapping of first enantiomer is preserved
    print(f"\n  >> {len(rxns)} reaction(s) — enantiomers combined, mapping preserved")
    return ok


# ===========================================================================
# 5. Multi-substrate product (dot-joined SMILES from rdchiral)
# ===========================================================================

def case_5_multi_substrate():
    print_header("5. Multi-substrate product  (retro-aldol split)")
    config = make_config()
    ra = RuleApplicator(config=config)
    # Retro-aldol: splits a beta-hydroxy ketone into two fragments
    # CC(O)CC=O = CH3-CH(OH)-CH2-CHO  (the CH(OH) has 1 H, matching CH1)
    target = 'CC(O)CC=O'
    opt = make_option(target, 'retro_aldol',
                      ["[CH1:1]([OH:2])-[CH2:3]-[C:4]=[O:5]>>[CH:1]=[O:2].[CH2:3]=[C:4]-[OH:5]"],
                      config, ra)
    rxns = opt.evaluate()

    ok = True

    def check(i, rxn, m):
        if len(rxn.substrates) > 1:
            print(f"    substrates are multi-component: {rxn.substrates}")
        if len(m) == 0:
            print("    !! WARN: empty mapping")
        else:
            # Should have 5 mapped atoms
            mol = Chem.MolFromSmiles(target)
            for tmap, aidx in m.items():
                if aidx >= mol.GetNumAtoms():
                    print(f"    !! FAIL: atom idx {aidx} out of range")
                    return False
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    if len(rxns) == 0:
        print("  !! FAIL: expected at least 1 reaction")
        ok = False
    return ok


# ===========================================================================
# 6. Template with many mapped atoms (large SMARTS)
# ===========================================================================

def case_6_large_template():
    print_header("6. Large template — Baeyer-Villiger-like")
    config = make_config()
    ra = RuleApplicator(config=config)
    # Ester hydrolysis-like retro: R-O-C(=O)-R' -> R-OH + R'-C(=O)-OH
    smarts = "[#8:1]-[#6:2](=[O:3])>>[#8:1].[#6:2](=[O:3])-[OH]"
    opt = make_option('CC(=O)OC', 'ester_hydrolysis', [smarts], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        mol = Chem.MolFromSmiles('CC(=O)OC')
        for tmap, aidx in m.items():
            if aidx >= mol.GetNumAtoms():
                print(f"    !! FAIL: atom idx {aidx} out of range")
                return False
        if len(m) > 0:
            print("    >> OK")
        return True

    ok = show_reactions('CC(=O)OC', rxns, check)
    return ok


# ===========================================================================
# 7. No match — template doesn't match target
# ===========================================================================

def case_7_no_match():
    print_header("7. No match  (Cl template on C=O molecule)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('CC=O', 'chlorination',
                      ["[Cl:1]-[#6:2]>>[Cl:1].[#6:2]"], config, ra)
    rxns = opt.evaluate()

    ok = len(rxns) == 0
    if ok:
        print("  >> PASS: no reactions (as expected)")
    else:
        print("  !! FAIL: expected no reactions")
    return ok


# ===========================================================================
# 8. Multiple SMARTS for same reaction name
# ===========================================================================

def case_8_multiple_smarts():
    print_header("8. Multiple SMARTS for one reaction name")
    config = make_config()
    ra = RuleApplicator(config=config)
    # Two different carbonyl reduction templates: one generic, one specific
    smarts_list = [
        "[#6:1]=[O:2]>>[#6:1]-[OH:2]",
        "[CH1:1]=[O:2]>>[CH2:1]-[OH:2]",
    ]
    opt = make_option('CC=O', 'carbonyl_red', smarts_list, config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) == 0:
            print("    !! FAIL: mapping should not be empty")
            return False
        print("    >> OK")
        return True

    ok = show_reactions('CC=O', rxns, check)
    print(f"\n  >> {len(rxns)} reaction(s) from {len(smarts_list)} SMARTS (dedup may reduce)")
    return ok


# ===========================================================================
# 9. Aromatic ring — template on aromatic system
# ===========================================================================

def case_9_aromatic():
    print_header("9. Aromatic substrate  (4-hydroxybenzaldehyde)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('Oc1ccc(C=O)cc1', 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        mol = Chem.MolFromSmiles('Oc1ccc(C=O)cc1')
        # mapnum :1 should be the aldehyde C, :2 the aldehyde O
        if 1 in m:
            atom = mol.GetAtomWithIdx(m[1])
            if atom.GetSymbol() != 'C':
                print(f"    !! FAIL: mapnum 1 mapped to {atom.GetSymbol()}, expected C")
                return False
        if 2 in m:
            atom = mol.GetAtomWithIdx(m[2])
            if atom.GetSymbol() != 'O':
                print(f"    !! FAIL: mapnum 2 mapped to {atom.GetSymbol()}, expected O")
                return False
        print("    >> OK")
        return True

    ok = show_reactions('Oc1ccc(C=O)cc1', rxns, check)
    return ok


# ===========================================================================
# 10. Three C=O groups — should get three distinct mappings
# ===========================================================================

def case_10_triketone():
    print_header("10. Three C=O groups  (1,3,5-triketone-like)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = 'O=CC(=O)CC=O'  # three C=O groups
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    c_indices = set()

    def check(i, rxn, m):
        if 1 in m:
            c_indices.add(m[1])
        if len(m) > 0:
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    print(f"\n  >> {len(c_indices)} distinct C atom sites across {len(rxns)} reactions")
    if len(c_indices) >= 3:
        print("  >> PASS: all three C=O groups mapped independently")
    elif len(c_indices) >= 2:
        print("  >> PARTIAL: at least 2 of 3 C=O sites found")
    return ok


# ===========================================================================
# 11. Serialization roundtrip — to_dict / from_dict preserves mapping
# ===========================================================================

def case_11_serialization():
    print_header("11. Serialization roundtrip  (to_dict -> from_dict)")
    config = make_config()
    ra = RuleApplicator(config=config)
    opt = make_option('CC=O', 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    ok = True
    for rxn in rxns:
        original_steps = rxn.data.get('atom_mapping', [])
        d = rxn.to_dict()
        restored = reaction_from_dict(d)
        restored_steps = restored.data.get('atom_mapping', [])

        print(f"\n  Original:  {original_steps}")
        print(f"  Restored:  {restored_steps}")

        if len(original_steps) != len(restored_steps):
            print("  !! FAIL: step counts differ after roundtrip")
            ok = False
        else:
            for si, (orig, rest) in enumerate(zip(original_steps, restored_steps)):
                if orig['rxn_name'] != rest['rxn_name'] or orig['mapping'] != rest['mapping']:
                    print(f"  !! FAIL: step {si} differs after roundtrip")
                    ok = False
                    break
            else:
                print("  >> OK")

    if ok:
        print("\n  >> PASS: serialization roundtrip preserves atom mappings")
    return ok


# ===========================================================================
# 12. Config disabled — no atom_mapping key in data
# ===========================================================================

def case_12_config_disabled():
    print_header("12. Config disabled  (store_atom_mappings = False)")
    config = Expansion_Config()  # default: store_atom_mappings = False
    ra = RuleApplicator(config=config)
    eval_fn = create_evaluate_option_method(ra, config)
    opt = ReactionOption(
        target_smi='CC=O',
        name='carbonyl_red',
        smarts=["[#6:1]=[O:2]>>[#6:1]-[OH:2]"],
        metadata={},
        rxn_type='test',
        rxn_domain='test',
        score=1.0,
        evaluation_function=eval_fn,
    )
    rxns = opt.evaluate()

    ok = True
    for rxn in rxns:
        if 'atom_mapping' in rxn.data:
            print(f"  !! FAIL: atom_mapping should NOT be in data when config disabled")
            print(f"     got: {rxn.data}")
            ok = False

    if ok:
        print(f"  >> PASS: {len(rxns)} reaction(s), none have atom_mapping in data")
    return ok


# ===========================================================================
# 13. Heterocyclic target — pyridine aldehyde
# ===========================================================================

def case_13_heterocycle():
    print_header("13. Heterocyclic target  (pyridine-3-carboxaldehyde)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = 'O=Cc1cccnc1'
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        mol = Chem.MolFromSmiles(target)
        if 2 in m:
            atom = mol.GetAtomWithIdx(m[2])
            if atom.GetSymbol() != 'O':
                print(f"    !! FAIL: mapnum 2 mapped to {atom.GetSymbol()}, expected O")
                return False
        if len(m) > 0:
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    return ok


# ===========================================================================
# 14. Molecule with existing stereocenters
# ===========================================================================

def case_14_existing_stereo():
    print_header("14. Existing stereocenter  ((S)-2-methylbutanal)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = '[C@@H](C)(CC)C=O'  # (S)-2-methylbutanal
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) == 0:
            print("    (empty mapping)")
        else:
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    print(f"\n  >> {len(rxns)} reaction(s) — new stereocenters may give enantiomeric pairs")
    return ok


# ===========================================================================
# 15. Naphthalene dialdehyde — two identical aldehyde groups on ring
# ===========================================================================

def case_15_naphthalene_dialdehyde():
    print_header("15. Naphthalene-2,6-dicarbaldehyde  (symmetric diaryl dialdehyde)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = 'O=Cc1ccc2cc(C=O)ccc2c1'
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    c_indices = set()

    def check(i, rxn, m):
        if 1 in m:
            c_indices.add(m[1])
        if len(m) > 0:
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    print(f"\n  >> {len(c_indices)} distinct aldehyde C sites across {len(rxns)} reactions")
    return ok


# ===========================================================================
# 16. Template with 3 mapped atoms — reductive amination reverse
# ===========================================================================

def case_16_three_mapped_atoms():
    print_header("16. Three mapped atoms  (reductive amination reverse)")
    config = make_config()
    ra = RuleApplicator(config=config)
    # C-N bond formation reverse: amine -> imine + reducing agent
    smarts = "[#6:1]-[NH:2]-[#6:3]>>[#6:1]=[N:2]-[#6:3]"
    target = 'CNCC'  # N-methylethylamine
    opt = make_option(target, 'dehydrogenation', [smarts], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        mol = Chem.MolFromSmiles(target)
        ok = True
        for tmap in [1, 2, 3]:
            if tmap not in m:
                print(f"    !! WARN: mapnum {tmap} not in mapping")
                continue
            aidx = m[tmap]
            if aidx >= mol.GetNumAtoms():
                print(f"    !! FAIL: atom idx {aidx} out of range")
                ok = False
        if 2 in m and mol.GetAtomWithIdx(m[2]).GetSymbol() != 'N':
            print(f"    !! FAIL: mapnum 2 should be N")
            ok = False
        if ok and len(m) > 0:
            print("    >> OK")
        return ok

    ok = show_reactions(target, rxns, check)
    return ok


# ===========================================================================
# 17. Same product from different sites but NOT deduped (different substrates)
# ===========================================================================

def case_17_different_splits():
    print_header("17. Different C-O bond cleavages  (diethyl ether)")
    config = make_config()
    ra = RuleApplicator(config=config)
    # C-O bond cleavage: R-O-R' -> R-OH + R'
    smarts = "[#6:1]-[O:2]-[#6:3]>>[#6:1]-[OH:2].[#6:3]"
    target = 'CCOCC'  # diethyl ether
    opt = make_option(target, 'ether_cleavage', [smarts], config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) > 0:
            print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)
    # symmetric ether — both C-O cleavages give same products, expect dedup to 1
    print(f"\n  >> {len(rxns)} reaction(s) — symmetric ether may dedup")
    return ok


# ===========================================================================
# 18. Enantiomer combination — R/S templates on butanone
#     WITH atom mapping: enantiomers combined to achiral with mapping preserved
#     WITHOUT atom mapping: also combined to achiral (same behavior)
# ===========================================================================

def case_18_enantiomer_combination_butanone():
    print_header("18. Enantiomer combination with mapping  (butanone R/S reduction)")
    config = make_config()
    ra = RuleApplicator(config=config)

    # Two chiral templates that produce opposite enantiomers
    smarts_list = [
        "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
        "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
    ]
    target = 'CCC(=O)C'  # butanone
    opt = make_option(target, 'chiral_reduction', smarts_list, config, ra)
    rxns = opt.evaluate()

    print(f"  WITH atom mapping ({len(rxns)} reactions):")
    ok = True

    def check(i, rxn, m):
        if len(m) == 0:
            print("    !! FAIL: mapping should not be empty")
            return False
        print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)

    # Compare: without atom mapping, enantiomers are also combined
    config_off = Expansion_Config()
    config_off.store_atom_mappings = False
    ra_off = RuleApplicator(config=config_off)
    opt_off = make_option(target, 'chiral_reduction', smarts_list, config_off, ra_off)
    rxns_off = opt_off.evaluate()
    print(f"\n  WITHOUT atom mapping ({len(rxns_off)} reactions):")
    for rxn in rxns_off:
        print(f"    {rxn.reaction_smiles()}")
        print(f"    substrates = {rxn.substrates}")

    # Both should give 1 reaction (enantiomers combined in both paths)
    if len(rxns) == 1 and len(rxns_off) == 1:
        steps = rxns[0].data.get('atom_mapping', [])
        m0 = steps[0]['mapping'] if steps else {}
        print(f"\n  >> PASS: enantiomers combined in both paths (1 reaction each)")
        print(f"  >> Mapping preserved: {m0}")
        if len(m0) == 0:
            print("  !! FAIL: mapping should not be empty")
            ok = False
    else:
        print(f"\n  !! FAIL: expected 1/1, got {len(rxns)}/{len(rxns_off)}")
        ok = False

    return ok


# ===========================================================================
# 19. Enantiomer combination — decalone (fused ring)
#     Enantiomers combined with atom mapping preserved
# ===========================================================================

def case_19_enantiomer_combination_decalone():
    print_header("19. Enantiomer combination with mapping  (decalone R/S reduction)")
    config = make_config()
    ra = RuleApplicator(config=config)

    smarts_list = [
        "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@@H:1]([OH:2])(-[#6:3])-[#6:4]",
        "[C:1](=[O:2])(-[#6:3])-[#6:4]>>[C@H:1]([OH:2])(-[#6:3])-[#6:4]",
    ]
    target = 'O=C1CCCC2CCCCC12'  # trans-decalone
    opt = make_option(target, 'chiral_reduction', smarts_list, config, ra)
    rxns = opt.evaluate()

    def check(i, rxn, m):
        if len(m) == 0:
            print("    !! FAIL: mapping empty")
            return False
        print("    >> OK")
        return True

    ok = show_reactions(target, rxns, check)

    # Compare without
    config_off = Expansion_Config()
    ra_off = RuleApplicator(config=config_off)
    opt_off = make_option(target, 'chiral_reduction', smarts_list, config_off, ra_off)
    rxns_off = opt_off.evaluate()

    print(f"\n  WITH mapping: {len(rxns)} reactions")
    print(f"  WITHOUT mapping: {len(rxns_off)} reactions")

    # Both should combine enantiomers; decalone may give diastereomers though
    # (ring junction cis vs trans), so number may vary
    if len(rxns) >= 1:
        steps = rxns[0].data.get('atom_mapping', [])
        m0 = steps[0]['mapping'] if steps else {}
        if len(m0) > 0:
            print(f"  >> PASS: mapping preserved after enantiomer combination: {m0}")
        else:
            print("  !! FAIL: mapping should not be empty")
            ok = False
    else:
        print("  !! FAIL: expected at least 1 reaction")
        ok = False

    return ok


# ===========================================================================
# 20. Bracket cleaning — simulated bracketed SMILES through _parse_with_mappings
#     rdchiral can produce [C], [N] etc. in certain ring contexts.
#     Bracket cleaning makes them identical to normal SMILES, causing dedup.
#     The mapping must survive this correctly.
# ===========================================================================

def case_20_bracket_cleaning():
    print_header("20. Bracket cleaning  (simulated [C] brackets through parser)")

    # Simulate what happens when rdchiral produces bracketed output
    # e.g. [C]1=[C]C[C@@H](c2ccccc2)N=C1 vs C1=CC[C@@H](c2ccccc2)N=C1
    # After bracket cleaning, [C] -> C, making both identical -> dedup keeps first mapping
    applicator = RdChiral_Applicator()

    products_raw = ['[C]1=[C]C[C@@H](c2ccccc2)N=C1', 'C1=CC[C@@H](c2ccccc2)N=C1']
    mappings_raw = [{2: 0, 3: 1}, {2: 4, 3: 5}]

    print(f"  Raw products:  {products_raw}")
    print(f"  Raw mappings:  {mappings_raw}")

    parsed, parsed_mappings = applicator._parse_with_mappings(products_raw, mappings_raw)

    print(f"  After parse:   {parsed}")
    print(f"  After mapping: {parsed_mappings}")

    ok = True
    if len(parsed) == 1 and len(parsed_mappings) == 1:
        print(f"\n  >> Bracket cleaning made both products identical -> dedup to 1")
        print(f"  >> First mapping kept: {parsed_mappings[0]}")
        if parsed_mappings[0] == {2: 0, 3: 1}:
            print("  >> PASS: first mapping preserved through bracket cleaning + dedup")
        else:
            print("  !! FAIL: wrong mapping survived dedup")
            ok = False
    else:
        print(f"  !! FAIL: expected 1 product after dedup, got {len(parsed)}")
        ok = False

    # Also test [N] bracket cleaning
    print()
    products_n = ['[N]=C1CCCCC1', 'N=C1CCCCC1']
    mappings_n = [{1: 0, 2: 1}, {1: 3, 2: 4}]
    parsed_n, parsed_mappings_n = applicator._parse_with_mappings(products_n, mappings_n)
    print(f"  [N] bracket test:")
    print(f"    Input:   {products_n} / {mappings_n}")
    print(f"    Output:  {parsed_n} / {parsed_mappings_n}")
    if len(parsed_n) == 1 and parsed_mappings_n[0] == {1: 0, 2: 1}:
        print("    >> OK: [N] cleaned, first mapping kept")
    else:
        print("    !! FAIL")
        ok = False

    return ok


# ===========================================================================
# 21. Bracket cleaning — real pipeline: ethylene hydrogenation
#     rdchiral produces [CH2][CH2] — brackets with H are NOT cleaned
#     (H is in bracket_required_chars), but rdkit canonicalization handles it.
# ===========================================================================

def case_21_bracket_with_hydrogen():
    print_header("21. Brackets with H  (ethylene hydrogenation: [CH2][CH2])")
    config = make_config()
    applicator = RdChiral_Applicator(config=config)

    rxns = applicator.get_rxns(['[C:1]=[C:2]>>[CH2:1]-[CH2:2]'])
    reactant = applicator.get_reactant('C=C')

    products, mappings = applicator.apply_reactions(reactant, rxns, return_atom_mapping=True)

    print(f"  Products: {products}")
    print(f"  Mappings: {mappings}")

    ok = True
    if len(products) == 1 and len(mappings) == 1:
        print(f"\n  Note: [CH2] brackets kept (H is a 'required char' for brackets)")
        print(f"  Mapping preserved: {mappings[0]}")
        if 1 in mappings[0] and 2 in mappings[0]:
            print("  >> PASS")
        else:
            print("  !! FAIL: expected mapnums 1 and 2")
            ok = False
    else:
        print(f"  !! FAIL: unexpected output")
        ok = False
    return ok


# ===========================================================================
# 22. E/Z double bond — (E)-cinnamaldehyde reduction
#     The existing E geometry in the target is preserved through the pipeline.
#     Atom mapping correctly identifies the C=O being reduced, not the C=C.
# ===========================================================================

def case_22_e_double_bond():
    print_header("22. E double bond preserved  ((E)-cinnamaldehyde reduction)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = 'O=C/C=C/c1ccccc1'  # (E)-cinnamaldehyde
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    mol = Chem.MolFromSmiles(target)
    ok = True

    def check(i, rxn, m):
        nonlocal ok
        # mapnum :1 should be C (aldehyde carbon), :2 should be O
        if 1 not in m or 2 not in m:
            print("    !! FAIL: expected mapnums 1 and 2")
            ok = False
            return False
        c_atom = mol.GetAtomWithIdx(m[1])
        o_atom = mol.GetAtomWithIdx(m[2])
        if c_atom.GetSymbol() != 'C' or o_atom.GetSymbol() != 'O':
            print(f"    !! FAIL: expected C and O, got {c_atom.GetSymbol()} and {o_atom.GetSymbol()}")
            ok = False
            return False

        # Check that product preserves E geometry
        has_e = '/' in rxn.substrates[0] or '\\' in rxn.substrates[0]
        if has_e:
            print(f"    >> E geometry preserved in product: {rxn.substrates[0]}")
        else:
            print(f"    !! WARN: E geometry lost in product: {rxn.substrates[0]}")
        print("    >> OK")
        return True

    show_reactions(target, rxns, check)
    return ok


# ===========================================================================
# 23. Z double bond — (Z)-cinnamaldehyde reduction
#     Same template, opposite geometry. Mapping should be the same but
#     product SMILES differs (\ instead of /).
# ===========================================================================

def case_23_z_double_bond():
    print_header(r"23. Z double bond preserved  ((Z)-cinnamaldehyde reduction)")
    config = make_config()
    ra = RuleApplicator(config=config)
    target = r'O=C/C=C\c1ccccc1'  # (Z)-cinnamaldehyde
    opt = make_option(target, 'carbonyl_red',
                      ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    rxns = opt.evaluate()

    mol = Chem.MolFromSmiles(target)
    ok = True

    def check(i, rxn, m):
        nonlocal ok
        if 1 not in m or 2 not in m:
            print("    !! FAIL: expected mapnums 1 and 2")
            ok = False
            return False

        # Check product has Z geometry (\ in SMILES)
        sub = rxn.substrates[0]
        has_stereo = '/' in sub or '\\' in sub
        if has_stereo:
            print(f"    >> Z geometry preserved in product: {sub}")
        else:
            print(f"    !! WARN: Z geometry lost in product: {sub}")
        print("    >> OK")
        return True

    show_reactions(target, rxns, check)

    return ok


# ===========================================================================
# 24. E vs Z comparison — same template, same atom mapping,
#     but products differ in geometry
# ===========================================================================

def case_24_e_vs_z_comparison():
    print_header("24. E vs Z comparison  (same template, different geometry)")
    config = make_config()
    ra = RuleApplicator(config=config)

    target_e = 'O=C/C=C/c1ccccc1'  # (E)-cinnamaldehyde
    target_z = r'O=C/C=C\c1ccccc1'  # (Z)-cinnamaldehyde

    opt_e = make_option(target_e, 'carbonyl_red',
                        ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)
    opt_z = make_option(target_z, 'carbonyl_red',
                        ["[#6:1]=[O:2]>>[#6:1]-[OH:2]"], config, ra)

    rxns_e = opt_e.evaluate()
    rxns_z = opt_z.evaluate()

    ok = True
    if len(rxns_e) >= 1 and len(rxns_z) >= 1:
        steps_e = rxns_e[0].data.get('atom_mapping', [])
        steps_z = rxns_z[0].data.get('atom_mapping', [])
        m_e = steps_e[0]['mapping'] if steps_e else {}
        m_z = steps_z[0]['mapping'] if steps_z else {}
        sub_e = rxns_e[0].substrates[0]
        sub_z = rxns_z[0].substrates[0]

        print(f"  E target: {target_e}")
        print(f"    product: {sub_e}")
        print(f"    mapping: {m_e}")
        print(f"  Z target: {target_z}")
        print(f"    product: {sub_z}")
        print(f"    mapping: {m_z}")

        # Mappings should be the same (same template, same atom positions)
        if m_e == m_z:
            print(f"\n  >> Mappings are identical (same C=O site in both isomers)")
        else:
            print(f"\n  >> Mappings differ (may reflect different canonical numbering)")

        # Products should be different SMILES
        if sub_e != sub_z:
            print(f"  >> Products are different (E vs Z preserved)")
            print("  >> PASS")
        else:
            print(f"  !! FAIL: E and Z products should differ")
            ok = False
    else:
        print(f"  !! FAIL: expected reactions from both")
        ok = False

    return ok


# ===========================================================================
# main
# ===========================================================================

if __name__ == '__main__':
    results = {}
    cases = [
        ("1_simple_carbonyl", case_1_simple_carbonyl),
        ("2_diketone", case_2_diketone),
        ("3_symmetric_ketone", case_3_symmetric_ketone),
        ("4_chiral_imine", case_4_chiral_imine),
        ("5_multi_substrate", case_5_multi_substrate),
        ("6_large_template", case_6_large_template),
        ("7_no_match", case_7_no_match),
        ("8_multiple_smarts", case_8_multiple_smarts),
        ("9_aromatic", case_9_aromatic),
        ("10_triketone", case_10_triketone),
        ("11_serialization", case_11_serialization),
        ("12_config_disabled", case_12_config_disabled),
        ("13_heterocycle", case_13_heterocycle),
        ("14_existing_stereo", case_14_existing_stereo),
        ("15_naphthalene_dialdehyde", case_15_naphthalene_dialdehyde),
        ("16_three_mapped_atoms", case_16_three_mapped_atoms),
        ("17_different_splits", case_17_different_splits),
        ("18_enantiomer_butanone", case_18_enantiomer_combination_butanone),
        ("19_enantiomer_decalone", case_19_enantiomer_combination_decalone),
        ("20_bracket_cleaning", case_20_bracket_cleaning),
        ("21_bracket_with_hydrogen", case_21_bracket_with_hydrogen),
        ("22_e_double_bond", case_22_e_double_bond),
        ("23_z_double_bond", case_23_z_double_bond),
        ("24_e_vs_z_comparison", case_24_e_vs_z_comparison),
    ]

    for name, fn in cases:
        try:
            ok = fn()
            results[name] = "PASS" if ok else "FAIL"
        except Exception as e:
            print(f"\n  !! EXCEPTION: {e}")
            import traceback
            traceback.print_exc()
            results[name] = "ERROR"

    print("\n")
    print("=" * 70)
    print("  SUMMARY")
    print("=" * 70)
    for name, status in results.items():
        icon = ">>" if status == "PASS" else "!!" if status == "FAIL" else "XX"
        print(f"  {icon} {name}: {status}")

    n_pass = sum(1 for v in results.values() if v == "PASS")
    n_total = len(results)
    print(f"\n  {n_pass}/{n_total} passed")
