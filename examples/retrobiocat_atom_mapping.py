"""
Example: RetroBioCat expander with atom mapping enabled.

Shows how to use config.store_atom_mappings = True with the RetroBioCat expander
to get step-based atom mapping on every Reaction object. Includes both
single-step and multi-step (two-step) reactions.

Run:  python examples/retrobiocat_atom_mapping.py
"""

from rdkit import Chem

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expanders.retrobiocat_expander import RetroBioCatExpander


def describe_step(step, indent=6):
    """Pretty-print one step of an atom mapping chain."""
    pad = " " * indent
    print(f"{pad}rxn_name:      {step['rxn_name']}")
    print(f"{pad}target_smi:    {step['target_smi']}")
    print(f"{pad}substrate_smis: {step['substrate_smis']}")
    mapping = step['mapping']
    print(f"{pad}mapping:       {mapping}")

    mol = Chem.MolFromSmiles(step['target_smi'])
    if mol and mapping:
        for tmap, aidx in sorted(mapping.items()):
            if aidx < mol.GetNumAtoms():
                sym = mol.GetAtomWithIdx(aidx).GetSymbol()
                print(f"{pad}  :{tmap} -> atom {aidx} ({sym})")


def example_single_step():
    """Single-step reactions: imine reduction on 2-phenylpiperidine."""
    print("=" * 70)
    print("  Single-step reactions: 2-phenylpiperidine imine")
    print("=" * 70)

    config = Expansion_Config()
    config.store_atom_mappings = True

    target = 'c1ccc(C2=NCCCC2)cc1'

    expander = RetroBioCatExpander(
        config=config,
        search_literature_precedent=False,
    )

    reactions = expander.get_reactions(target)

    print(f"\nTarget: {target}")
    print(f"Reactions found: {len(reactions)}\n")

    for i, rxn in enumerate(reactions):
        print(f"  Reaction {i+1}: {rxn.name}")
        print(f"    {rxn.product} -> {rxn.substrates}")
        steps = rxn.data.get('atom_mapping', [])
        print(f"    atom_mapping: {len(steps)} step(s)")
        for si, step in enumerate(steps):
            print(f"    --- step {si+1} ---")
            describe_step(step)
        print()


def example_multi_step():
    """Multi-step (two-step) reactions: secondary amine reductions etc."""
    print("=" * 70)
    print("  Multi-step reactions: secondary amine target")
    print("=" * 70)

    config = Expansion_Config()
    config.store_atom_mappings = True

    # This target can trigger two-step reactions like "Secondary amine reduction"
    # which is "Imine reduction" then "Amine oxidation"
    target = 'c1ccc([C@@H]2CCCCN2)cc1'

    expander = RetroBioCatExpander(
        config=config,
        search_literature_precedent=False,
        include_two_step=True,
    )

    reactions = expander.get_reactions(target)

    print(f"\nTarget: {target}")
    print(f"Reactions found: {len(reactions)}\n")

    # Separate single-step and multi-step reactions
    single_step = [r for r in reactions if len(r.data.get('atom_mapping', [])) == 1]
    multi_step = [r for r in reactions if len(r.data.get('atom_mapping', [])) > 1]
    no_mapping = [r for r in reactions if 'atom_mapping' not in r.data]

    print(f"  Single-step with mapping: {len(single_step)}")
    print(f"  Multi-step with mapping:  {len(multi_step)}")
    print(f"  Without mapping:          {len(no_mapping)}")
    print()

    # Show single-step examples (first 3)
    if single_step:
        print("  --- Single-step examples ---")
        for rxn in single_step[:3]:
            step = rxn.data['atom_mapping'][0]
            print(f"    {rxn.name}: {rxn.product} -> {rxn.substrates}")
            print(f"      mapping: {step['mapping']}")
        print()

    # Show all multi-step reactions
    if multi_step:
        print("  --- Multi-step reactions ---")
        for rxn in multi_step:
            steps = rxn.data['atom_mapping']
            print(f"\n    {rxn.name} ({len(steps)} steps)")
            print(f"    overall: {rxn.product} -> {rxn.substrates}")
            for si, step in enumerate(steps):
                print(f"    step {si+1}:")
                describe_step(step, indent=6)

            # Verify chain property: step[0].target_smi == original target
            assert steps[0]['target_smi'] == rxn.product, "First step target should be original target"
            # Verify chain property: step[-1].substrate_smis == final substrates
            assert steps[-1]['substrate_smis'] == rxn.substrates, "Last step substrates should be final substrates"
            print(f"    chain properties verified")
    else:
        print("  (no multi-step reactions found for this target)")
    print()


def example_config_off():
    """Confirm that other expanders (or config off) produce no atom_mapping."""
    print("=" * 70)
    print("  Config disabled: no atom_mapping in data")
    print("=" * 70)

    config = Expansion_Config()
    # config.store_atom_mappings defaults to False

    target = 'c1ccc(C2=NCCCC2)cc1'

    expander = RetroBioCatExpander(
        config=config,
        search_literature_precedent=False,
    )

    reactions = expander.get_reactions(target)

    print(f"\nTarget: {target}")
    print(f"config.store_atom_mappings = {config.store_atom_mappings}")
    print(f"Reactions found: {len(reactions)}")

    has_mapping = sum(1 for r in reactions if 'atom_mapping' in r.data)
    print(f"Reactions with atom_mapping: {has_mapping}")
    assert has_mapping == 0, "No reactions should have atom_mapping when config is off"
    print("PASS: no atom_mapping keys present\n")


if __name__ == '__main__':
    example_single_step()
    example_multi_step()
    example_config_off()
