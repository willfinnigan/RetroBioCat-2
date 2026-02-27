# Atom Mapping

Atom mapping tracks which atoms in the **target molecule** correspond to which
atoms in the **reaction template**. When enabled, every `Reaction` object
carries a step-based mapping in `reaction.data['atom_mapping']`.

## Enabling atom mapping

```python
from rbc2.configs.expansion_config import Expansion_Config
from rbc2.expansion.expanders.retrobiocat_expander import RetroBioCatExpander

config = Expansion_Config()
config.store_atom_mappings = True

expander = RetroBioCatExpander(config=config)
reactions = expander.get_reactions('c1ccc(C2=NCCCC2)cc1')

for rxn in reactions:
    for step in rxn.data['atom_mapping']:
        print(step['rxn_name'], step['mapping'])
```

Set `config.store_atom_mappings = True` to enable. Only the `RetroBioCatExpander`
supports atom mapping — other expanders ignore this setting.

When using `MultiExpander`, all expanders share one config. Setting the flag
only activates mapping for RetroBioCat; the others are unaffected.

## Data format

`reaction.data['atom_mapping']` is always a **list of step dicts**, even for
single-step reactions. This gives a consistent interface for both single-step
and multi-step reactions.

Each step dict has four keys:

| Key | Type | Description |
|-----|------|-------------|
| `rxn_name` | `str` | Human-readable reaction name (e.g. "Imine reduction") |
| `target_smi` | `str` | The molecule this step acts on |
| `substrate_smis` | `List[str]` | The substrates produced by this step |
| `mapping` | `dict` | `{template_mapnum: target_atom_idx}` — see below |

### Single-step example

```python
reaction.data['atom_mapping'] = [
    {
        'rxn_name': 'Ketone reduction',
        'target_smi': 'CC=O',
        'substrate_smis': ['CCO'],
        'mapping': {1: 1, 2: 2}
    }
]
```

### Multi-step example

Multi-step reactions (e.g. "Amine deracemization" = "Imine reduction" then
"Amine oxidation") produce a chain of step dicts:

```python
reaction.data['atom_mapping'] = [
    {
        'rxn_name': 'Imine reduction',
        'target_smi': 'c1ccc([C@@H]2CCCCN2)cc1',
        'substrate_smis': ['c1ccc(C2=NCCCC2)cc1'],
        'mapping': {1: 5, 3: 4, 4: 9}
    },
    {
        'rxn_name': 'Amine oxidation',
        'target_smi': 'c1ccc(C2=NCCCC2)cc1',
        'substrate_smis': ['c1ccc(C2CCCCN2)cc1'],
        'mapping': {1: 4, 2: 5, 3: 6}
    }
]
```

**Chain properties:**
- `steps[0]['target_smi']` is the original target molecule
- `steps[-1]['substrate_smis']` is the final substrates (= `reaction.substrates`)
- Each step's `substrate_smis` feeds into the next step's `target_smi`

## Reading the mapping

Each `mapping` dict maps **template atom map number → target molecule atom
index (0-based)**. Template map numbers come from the SMARTS pattern — e.g.
in `[C:1]=[O:2]>>[C:1][OH:2]`, `:1` is the carbon and `:2` is the oxygen.

You can resolve these to atom symbols with RDKit:

```python
from rdkit import Chem

step = reaction.data['atom_mapping'][0]
mol = Chem.MolFromSmiles(step['target_smi'])

for mapnum, atom_idx in step['mapping'].items():
    atom = mol.GetAtomWithIdx(atom_idx)
    print(f"template :{mapnum} -> atom {atom_idx} ({atom.GetSymbol()})")
```

## Configuration

Two things must be true for atom mapping to be produced:

| Setting | Where | Default | Purpose |
|---------|-------|---------|---------|
| `config.store_atom_mappings` | `Expansion_Config` | `False` | User opt-in for atom mapping |
| `expander.supports_atom_mapping` | Expander class | `False` | Expander declares capability |

Currently only `RetroBioCatExpander` sets `supports_atom_mapping = True`.
Other expanders (AIZynthfinder, RingBreaker, etc.) leave it `False`.

## Data flow

```
Expansion_Config.store_atom_mappings = True
                    +
Expander.supports_atom_mapping = True
                    │
                    ▼
         RetroBioCatExpander._get_reactions()
                    │
                    ▼
         RuleApplicator.apply(smi, rxns, return_atom_mapping=True)
            │                               │
            ▼                               ▼
   RdChiral_Applicator              Multi_step_applicator
   .apply_rules()                   .apply_multi_step_rules()
            │                               │
            ▼                               ▼
      rdchiralRun()                 _apply_steps_with_mappings()
   {smi: {mapnum: idx}}            builds chain of step dicts
            │                               │
            └──────────┬────────────────────┘
                       ▼
              create_reactions()
                       │
            ┌──────────┴──────────┐
            ▼                     ▼
     Single-step:            Multi-step:
     wraps flat dict         passes chain
     in 1-element list       through as-is
            │                     │
            └──────────┬──────────┘
                       ▼
         Reaction.data['atom_mapping']
              List[step_dict]
```

## Key files

| File | Role |
|------|------|
| `rbc2/template_application/apply_template/rdchiral/main.py` | `rdchiralRun` — extracts `{template_mapnum: mol_atom_idx}` from rdchiral's internal atom tracking |
| `rbc2/template_application/apply_template/components/rdchiral_applicator.py` | Runs rdchiral for each template, carries mappings through parsing/dedup/enantiomer-combination |
| `rbc2/template_application/apply_template/components/multi_step_applicator.py` | `_apply_steps_with_mappings()` builds the chain of step dicts for multi-step reactions |
| `rbc2/template_application/apply_template/rule_applicator.py` | Orchestrates single-step and multi-step, threads `return_atom_mapping` to both |
| `rbc2/template_application/create_reactions_from_output/create_reactions.py` | Wraps single-step mappings in step format, passes multi-step chains through |
| `rbc2/expansion/expanders/retrobiocat_expander.py` | Declares `supports_atom_mapping = True`, calls the pipeline |
| `rbc2/expansion/default_expander.py` | Declares `supports_atom_mapping = False` (base class default) |
| `rbc2/configs/expansion_config.py` | `store_atom_mappings` flag |

## Enantiomer handling

When `config.combine_enantiomers = True` (default), R/S enantiomers produced
by chiral templates are merged into a single achiral product. The mapping from
the first enantiomer is preserved. This is safe because enantiomers from the
same template match always map to the same target atoms.

## When atom mapping is absent

- `config.store_atom_mappings = False` (the default) — `reaction.data` has no
  `'atom_mapping'` key
- Expanders other than RetroBioCat — same, no key present
- Use `rxn.data.get('atom_mapping', [])` to safely handle both cases

