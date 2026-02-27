"""Augmented Transformer model for template-free retrosynthesis.

Wraps OpenNMT-py's ServerModel for inference. OpenNMT-py is an optional
dependency - the model will raise ImportError if it's not installed.
"""
import glob
import logging
import os
import re
from typing import List, Optional

import numpy as np
from rdkit import Chem

from rbc2.configs.data_path import path_to_data_folder
from rbc2.expansion.expanders.template_free.models.model_interface import (
    TemplateFreeModel, TemplateFreeResult
)

logger = logging.getLogger(__name__)


def smi_tokenizer(smi: str) -> str:
    """Tokenize a SMILES string into space-separated tokens."""
    pattern = (r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|"
               r"\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])")
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    return " ".join(tokens)


def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize a SMILES string, removing atom mapping."""
    smiles = "".join(smiles.split())
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    for a in mol.GetAtoms():
        if a.HasProp('molAtomMapNumber'):
            a.ClearProp('molAtomMapNumber')
    cano = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    # Double canonicalization for stereochemistry edge cases
    mol2 = Chem.MolFromSmiles(cano)
    if mol2 is not None:
        cano = Chem.MolToSmiles(mol2, isomericSmiles=True, canonical=True)
    return cano


class AugmentedTransformerModel(TemplateFreeModel):
    """Inference wrapper for the Augmented Transformer retrosynthesis model.

    Uses OpenNMT-py's ServerModel for beam search decoding.
    Requires: pip install OpenNMT-py
    """

    def __init__(self,
                 model_dir: Optional[str] = None,
                 beam_size: int = 30,
                 n_best: int = 30):
        if model_dir is None:
            model_dir = os.path.join(path_to_data_folder, 'template_free',
                                     'augmented_transformer')
        self.model_dir = model_dir
        self.beam_size = beam_size
        self.n_best = n_best
        self._model = None
        self._loaded = False

    def is_loaded(self) -> bool:
        return self._loaded

    def load(self):
        try:
            from onmt.translate.translation_server import ServerModel as ONMTServerModel
        except ImportError:
            raise ImportError(
                "OpenNMT-py is required for AugmentedTransformerModel. "
                "Install it with: pip install 'OpenNMT-py==2.3.0'"
            )

        from rbc2.configs.download_data_files.download_template_free import ensure_augmented_transformer_exists
        ensure_augmented_transformer_exists()

        # Find checkpoint
        checkpoint_file = os.path.join(self.model_dir, "model_step_500000.pt")
        if not os.path.isfile(checkpoint_file):
            checkpoint_list = sorted(glob.glob(
                os.path.join(self.model_dir, "model_step_*.pt")))
            if not checkpoint_list:
                raise FileNotFoundError(
                    f"No model checkpoint found in {self.model_dir}. "
                    f"Place model_step_*.pt file(s) in this directory."
                )
            checkpoint_file = checkpoint_list[-1]
            logger.info("Using checkpoint: %s", checkpoint_file)

        import torch
        use_cpu = not torch.cuda.is_available()

        onmt_config = {
            "models": checkpoint_file,
            "n_best": self.n_best,
            "beam_size": self.beam_size,
            "gpu": -1 if use_cpu else 0
        }

        self._model = ONMTServerModel(
            opt=onmt_config,
            model_id=0,
            load=True
        )
        self._loaded = True
        logger.info("Augmented Transformer model loaded from %s", checkpoint_file)

    def predict(self, smi: str, num_results: int = 10) -> List[TemplateFreeResult]:
        if not self.is_loaded():
            self.load()

        # Canonicalize and tokenize input
        canon_smi = canonicalize_smiles(smi)
        if not canon_smi:
            return []

        tokenized = smi_tokenizer(canon_smi)
        inputs = [{"src": tokenized}]

        # Run inference
        try:
            products, scores, _, _, _ = self._model.run(inputs=inputs)
        except Exception as e:
            logger.warning("Augmented Transformer inference failed for %s: %s",
                           smi, e)
            return []

        # Post-process: deduplicate, canonicalize, aggregate scores
        valid_product_scores = {}
        for prod, score in zip(products, scores):
            product = "".join(prod.split())
            product = canonicalize_smiles(product)
            if not product:
                continue

            # Filter predictions that lose too many atoms (>4)
            try:
                input_mol = Chem.MolFromSmiles(canon_smi)
                prod_mol = Chem.MolFromSmiles(product)
                if input_mol and prod_mol:
                    if (input_mol.GetNumAtoms() - prod_mol.GetNumAtoms()) > 4:
                        continue
            except Exception:
                continue

            # Aggregate scores for duplicate predictions
            exp_score = np.exp(score)
            if product in valid_product_scores:
                valid_product_scores[product] += exp_score
            else:
                valid_product_scores[product] = exp_score

        # Sort by score descending
        sorted_products = sorted(valid_product_scores.items(),
                                 key=lambda x: x[1], reverse=True)

        # Normalize scores by total probability mass to get [0, 1] range
        total_score = sum(valid_product_scores.values())
        if total_score == 0:
            total_score = 1.0

        results = []
        for reactants_smi, score in sorted_products[:num_results]:
            norm_score = score / total_score
            results.append(TemplateFreeResult(
                reactants=reactants_smi,
                score=norm_score
            ))

        return results
