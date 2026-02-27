"""Base class for graph-based template-free inference models.

Extracts shared logic from EGraph2RSmilesModel and Graph2SmilesModel:
vocabulary loading, checkpoint loading, beam search decoding, and token
decoding to SMILES strings.
"""

import logging
import math
import os
from abc import abstractmethod
from typing import List, Optional

import torch

from rbc2.expansion.expanders.template_free.models.model_interface import TemplateFreeModel, TemplateFreeResult
from rbc2.expansion.expanders.template_free.models.shared.beam_search import beam_search_decode
from rbc2.expansion.expanders.template_free.utils.data_utils import load_vocab_from_file

logger = logging.getLogger(__name__)


def get_best_device() -> torch.device:
    """Select the best device for inference.

    Defaults to CPU. Single-molecule beam search on these small models
    (256 hidden dim) is faster on CPU than on GPU (CUDA or MPS) because
    the sequential decoding steps and CPU-GPU transfer overhead outweigh
    any parallelism gains. GPU would only help with batched inference or
    larger models.
    """
    return torch.device("cpu")


class GraphInferenceModel(TemplateFreeModel):
    """Base inference wrapper for graph-encoder + transformer-decoder models.

    Subclasses must implement:
    - _ensure_weights_exist(): Download model weights if missing.
    - _create_config(): Return a default config dataclass instance.
    - _create_model(config, vocab): Instantiate the nn.Module.
    - _create_batch(smi): Build input batch from SMILES.
    - _normalize_score(beam_score): Convert beam search log-prob to [0, 1].
    """

    def __init__(self,
                 checkpoint_path: str,
                 vocab_path: str,
                 beam_size: int = 10,
                 max_length: int = 512,
                 temperature: float = 1.0):
        self.checkpoint_path = checkpoint_path
        self.vocab_path = vocab_path
        self.beam_size = beam_size
        self.max_length = max_length
        self.temperature = temperature

        self._model = None
        self._vocab = None
        self._vocab_tokens = None
        self._device = None

    def is_loaded(self) -> bool:
        return self._model is not None

    def load(self):
        """Load model weights and vocabulary."""
        self._ensure_weights_exist()

        self._device = get_best_device()

        self._vocab = load_vocab_from_file(self.vocab_path)
        self._vocab_tokens = {v: k for k, v in self._vocab.items()}

        config = self._create_config()
        state = torch.load(self.checkpoint_path, map_location="cpu", weights_only=False)

        if "args" in state:
            checkpoint_args = state["args"]
            for k, v in vars(checkpoint_args).items():
                if hasattr(config, k):
                    setattr(config, k, v)
            self._on_config_loaded(config)

        self._model = self._create_model(config, self._vocab)

        state_dict = state.get("state_dict", state)
        state_dict = {k.replace("module.", ""): v for k, v in state_dict.items()}
        self._model.load_state_dict(state_dict, strict=False)
        self._model.to(self._device)
        self._model.eval()

        logger.info(f"{self.__class__.__name__} loaded on {self._device}")

    def predict(self, smi: str, num_results: int = 10) -> List[TemplateFreeResult]:
        """Predict reactants for a product SMILES."""
        if not self.is_loaded():
            self.load()

        batch = self._create_batch(smi)

        with torch.no_grad():
            results = beam_search_decode(
                model=self._model,
                batch=batch,
                vocab=self._vocab,
                beam_size=self.beam_size,
                n_best=num_results,
                max_length=self.max_length,
                temperature=self.temperature,
                device=self._device)

        template_free_results = []
        for result in results:
            tokens = result.token_ids
            if tokens and tokens[0] == self._vocab.get("_SOS", 0):
                tokens = tokens[1:]
            if tokens and tokens[-1] == self._vocab.get("_EOS", 0):
                tokens = tokens[:-1]

            smi_tokens = []
            for tid in tokens:
                token = self._vocab_tokens.get(tid, "")
                if token and not token.startswith("_"):
                    smi_tokens.append(token)

            reactants_smi = "".join(smi_tokens)
            if reactants_smi:
                template_free_results.append(TemplateFreeResult(
                    reactants=reactants_smi,
                    score=self._normalize_score(result.score)
                ))

        return template_free_results

    @abstractmethod
    def _ensure_weights_exist(self):
        """Download model weights if not present."""
        pass

    @abstractmethod
    def _create_config(self):
        """Return a default config dataclass instance."""
        pass

    @abstractmethod
    def _create_model(self, config, vocab: dict) -> torch.nn.Module:
        """Instantiate the neural network model."""
        pass

    @abstractmethod
    def _create_batch(self, smi: str):
        """Create input batch from a SMILES string."""
        pass

    def _normalize_score(self, beam_score: float) -> float:
        """Convert beam search score (length-normalized log-prob) to [0, 1].

        Default: exp(log_prob) which maps log-probs to probabilities.
        Subclasses may override for model-specific normalization.
        """
        return math.exp(beam_score)

    def _on_config_loaded(self, config):
        """Hook called after loading config from checkpoint. Override if needed."""
        pass
