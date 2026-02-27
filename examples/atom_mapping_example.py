"""
Example: Exposing atom mapping from rdchiral template application.

This script demonstrates how to use the return_atom_mapping feature at four levels:
1. Low-level: rdchiralRun directly
2. Mid-level: RdChiral_Applicator
3. High-level: RuleApplicator
4. Full pipeline: ReactionOption.evaluate() with atom mappings stored in Reaction.data
"""

from rdkit import Chem

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.reaction_option import ReactionOption, create_evaluate_option_method
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction, rdchiralReactants
from rbc2.template_application.apply_template.rdchiral.main import rdchiralRun
from rbc2.template_application.apply_template.components.rdchiral_applicator import RdChiral_Applicator
from rbc2.template_application.apply_template.rule_applicator import RuleApplicator


def example_1_rdchiral_run():
    """Low-level: rdchiralRun with atom mapping on a simple carbonyl reduction."""
    print("=" * 60)
    print("Example 1: rdchiralRun — carbonyl reduction of acetaldehyde")
    print("=" * 60)

    smarts = '[C:1]=[O:2]>>[C:1][OH:2]'
    target_smi = 'CC=O'  # acetaldehyde

    rxn = rdchiralReaction(smarts)
    reactants = rdchiralReactants(target_smi)

    outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

    print(f"Target molecule: {target_smi}")
    print(f"Template SMARTS: {smarts}")
    print(f"Outcomes: {outcomes}")
    print()

    # Show atom labels for reference
    mol = Chem.MolFromSmiles(target_smi)
    print("Target molecule atoms:")
    for atom in mol.GetAtoms():
        print(f"  idx={atom.GetIdx()}: {atom.GetSymbol()} (in SMILES: {target_smi})")

    print()
    for smi, mapping in atom_mappings.items():
        print(f"Outcome: {smi}")
        print(f"  Atom mapping (template_mapnum -> target_atom_idx):")
        for template_num, atom_idx in sorted(mapping.items()):
            atom = mol.GetAtomWithIdx(atom_idx)
            print(f"    template :{template_num} -> atom idx {atom_idx} ({atom.GetSymbol()})")
    print()


def example_2_diketone():
    """Low-level: multiple matches on a diketone — each C=O maps to different atoms."""
    print("=" * 60)
    print("Example 2: rdchiralRun — diketone with two C=O groups")
    print("=" * 60)

    smarts = '[C:1]=[O:2]>>[C:1][OH:2]'
    target_smi = 'O=CC(=O)C'  # methylglyoxal

    rxn = rdchiralReaction(smarts)
    reactants = rdchiralReactants(target_smi)

    outcomes, atom_mappings = rdchiralRun(rxn, reactants, return_atom_mapping=True)

    mol = Chem.MolFromSmiles(target_smi)

    print(f"Target molecule: {target_smi}")
    print(f"Template SMARTS: {smarts}")
    print(f"Number of outcomes: {len(outcomes)}")
    print()

    print("Target molecule atoms:")
    for atom in mol.GetAtoms():
        print(f"  idx={atom.GetIdx()}: {atom.GetSymbol()}")
    print()

    for i, smi in enumerate(outcomes):
        mapping = atom_mappings[smi]
        print(f"Outcome {i+1}: {smi}")
        print(f"  Mapping: {mapping}")
        c_idx = mapping[1]  # template :1 is the C
        o_idx = mapping[2]  # template :2 is the O
        print(f"  C atom (template :1) -> target atom idx {c_idx} ({mol.GetAtomWithIdx(c_idx).GetSymbol()})")
        print(f"  O atom (template :2) -> target atom idx {o_idx} ({mol.GetAtomWithIdx(o_idx).GetSymbol()})")
    print()


def example_3_applicator():
    """Mid-level: using RdChiral_Applicator with atom mapping (parallel lists)."""
    print("=" * 60)
    print("Example 3: RdChiral_Applicator — multiple reaction rules")
    print("=" * 60)

    applicator = RdChiral_Applicator()

    rxns_dict = {
        'carbonyl_reduction': applicator.get_rxns(["[#6:1]=[O:2]>>[#6:1]-[OH:2]"]),
        'imine_reduction': applicator.get_rxns(["[#6:1]=[N:2]>>[#6:1]-[NH:2]"]),
    }

    target_smi = 'O=CC=N'  # has both C=O and C=N

    products, atom_mappings = applicator.apply_rules(target_smi, rxns_dict,
                                                      return_atom_mapping=True)

    print(f"Target molecule: {target_smi}")
    print(f"Rules applied: {list(rxns_dict.keys())}")
    print()

    mol = Chem.MolFromSmiles(target_smi)
    for rxn_name, product_list in products.items():
        mapping_list = atom_mappings.get(rxn_name, [])
        print(f"Reaction: {rxn_name}")
        for i, product_set in enumerate(product_list):
            mapping = mapping_list[i] if i < len(mapping_list) else {}
            print(f"  Products: {product_set}")
            print(f"  Mapping: {mapping}")
    print()


def example_4_rule_applicator():
    """High-level: using RuleApplicator with atom mapping (parallel dicts)."""
    print("=" * 60)
    print("Example 4: RuleApplicator — full pipeline with atom mapping")
    print("=" * 60)

    applicator = RuleApplicator()

    rxns = {
        'ketone_reduction': applicator.smarts_to_rdchiral(
            ["[#6:1](-[#6:2])=[O:3]>>[#6:1](-[#6:2])-[OH:3]"]
        ),
    }

    target_smi = 'CC(=O)CC(=O)C'  # 2,4-pentanedione

    products, atom_mappings = applicator.apply(target_smi, rxns, return_atom_mapping=True)

    mol = Chem.MolFromSmiles(target_smi)

    print(f"Target molecule: {target_smi}")
    print(f"Products: {products}")
    print()

    print("Target molecule atoms:")
    for atom in mol.GetAtoms():
        print(f"  idx={atom.GetIdx()}: {atom.GetSymbol()}")
    print()

    print("Atom mappings per outcome (parallel lists):")
    for rxn_name, product_list in products.items():
        mapping_list = atom_mappings.get(rxn_name, [])
        for i, product_set in enumerate(product_list):
            mapping = mapping_list[i] if i < len(mapping_list) else {}
            print(f"  Outcome: {product_set}")
            for tmap, aidx in sorted(mapping.items()):
                atom = mol.GetAtomWithIdx(aidx)
                print(f"    template :{tmap} -> atom idx {aidx} ({atom.GetSymbol()})")
    print()


def example_5_reaction_option():
    """Full pipeline: ReactionOption.evaluate() stores atom mappings in Reaction.data."""
    print("=" * 60)
    print("Example 5: ReactionOption — atom mappings on Reaction objects")
    print("=" * 60)

    config = Expansion_Config()
    config.store_atom_mappings = True

    rule_applicator = RuleApplicator(config=config)
    eval_fn = create_evaluate_option_method(rule_applicator, config)

    option = ReactionOption(
        target_smi='O=CC(=O)C',  # methylglyoxal — two C=O groups
        name='carbonyl_reduction',
        smarts=["[#6:1]=[O:2]>>[#6:1]-[OH:2]"],
        metadata={},
        rxn_type='test',
        rxn_domain='test',
        score=1.0,
        evaluation_function=eval_fn,
    )

    reactions = option.evaluate()

    mol = Chem.MolFromSmiles(option.target_smi)

    print(f"Target molecule: {option.target_smi}")
    print(f"Number of reactions: {len(reactions)}")
    print()

    for i, rxn in enumerate(reactions):
        steps = rxn.data.get('atom_mapping', [])
        print(f"Reaction {i+1}: {rxn.reaction_smiles()}")
        print(f"  atom_mapping (step-based): {len(steps)} step(s)")
        for step_idx, step in enumerate(steps):
            print(f"  Step {step_idx+1}:")
            print(f"    rxn_name: {step['rxn_name']}")
            print(f"    target_smi: {step['target_smi']}")
            print(f"    substrate_smis: {step['substrate_smis']}")
            print(f"    mapping: {step['mapping']}")
            step_mol = Chem.MolFromSmiles(step['target_smi'])
            if step_mol and step['mapping']:
                for tmap, aidx in sorted(step['mapping'].items()):
                    if aidx < step_mol.GetNumAtoms():
                        atom = step_mol.GetAtomWithIdx(aidx)
                        print(f"      template :{tmap} -> atom idx {aidx} ({atom.GetSymbol()})")
    print()


if __name__ == '__main__':
    example_1_rdchiral_run()
    example_2_diketone()
    example_3_applicator()
    example_4_rule_applicator()
    example_5_reaction_option()
