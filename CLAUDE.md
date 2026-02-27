# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RetroBioCat 2 (`rbc2`) is a Python package for computer-aided retrosynthesis planning combining biocatalysis and chemistry. It uses Monte Carlo Tree Search (MCTS) to find multi-step synthesis routes, with pluggable "expander" backends for different reaction domains.

## Build & Development Commands

```bash
# Install in editable mode
pip install -e .

# Run all tests
pytest tests

# Run a single test file
pytest tests/test_expansion/test_retrobiocat_expander.py

# Run a specific test
pytest tests/test_expansion/test_retrobiocat_expander.py::test_function_name

# Run tests with coverage (as CI does)
pytest tests --cov=rbc2 --cov-report=html
```

**Python version:** 3.9+

**ARM Mac note:** HDF5/tables requires homebrew exports before pip install:
```bash
export HDF5_DIR=$(brew --prefix hdf5)
export BLOSC_DIR=$(brew --prefix c-blosc)
pip install tables
```

## Architecture

### Data Model (`rbc2/data_model/`)
- **Reaction**: Core dataclass — product SMILES, substrate list, score, rxn_type, rxn_domain, metadata, precedents
- **ReactionOption**: Unexpanded template (lazy evaluation). Calling `.evaluate()` applies SMARTS templates to produce Reactions
- **Network**: Graph tracking all discovered molecules and reactions. Prevents re-expansion of already-seen molecules
- **Pathway**: Linear sequence of Reactions representing a synthesis route

### Expander System (`rbc2/expansion/`)
Plugin architecture for retrosynthesis backends. All expanders implement the `Expander` ABC (`expander_interface.py`):
- `get_reactions(smi)` → `List[Reaction]`
- `get_options(smi)` → `List[ReactionOption]`

Available expanders (registered in `expander_repository.py`):
- **retrobiocat** — curated enzyme reaction rules (biocatalysis)
- **aizynthfinder**, **ringbreaker**, **askcos** — neural network policy-based (chemistry)
- **retrorules**, **bkms**, **enzymemap** — template-based (biosynthesis)

The `DefaultExpander` base class handles the common pipeline: create options → apply templates via `RuleApplicator` → process reactions → search precedents.

`MultiExpander` combines multiple expanders for MCTS, supporting interleaving or score-based ordering.

### MCTS (`rbc2/mcts/`)
Core search algorithm in `mcts.py`. Loop: **selection** (UCB1) → **expansion** (apply expanders) → **rollout** (depth-first) → **scoring** → **backpropagation**. Solved pathways are those where all leaf molecules are buyable starting materials.

### Template Application (`rbc2/template_application/`)
Applies SMARTS reaction templates to molecules. Two backends: RDKit (default off) and RDChiral (default, handles stereochemistry). Post-processing includes deduplication and sanity checking.

### Reaction Evaluation (`rbc2/reaction_evaluation/`)
- **Feasibility filters**: AIZynthfinder neural network filter, RetroBioCat precedent-based filter
- **Starting material evaluators**: `CommercialSME` (commercial availability), `EcoliSME` (E. coli metabolism)
- **Complexity scoring**: SCScore model for synthetic complexity

### Configuration (`rbc2/configs/`)
Two main config classes, both supporting `update_from_dict()` / `to_dict()` serialization:
- **Expansion_Config**: Rule application settings (chiral handling, deduplication, rdchiral vs rdkit)
- **MCTS_Config**: Search parameters (exploration constant, time/iteration limits, scoring mode, domain restrictions)

Data file download scripts live in `rbc2/configs/download_data_files/`.

## Key Patterns

- **SMILES everywhere**: Molecules are represented as SMILES strings (variable name typically `smi`)
- **Lazy evaluation**: ReactionOption holds template info; Reactions are only created when `.evaluate()` is called
- **Config objects**: Passed through the call chain; created with defaults, overridden via constructor or `update_from_dict()`
- **Network caching**: The Network object is shared across expanders to avoid redundant computation
- **Expander plugin registration**: `get_expanders(['retrobiocat', 'aizynthfinder'])` returns instantiated expander list
