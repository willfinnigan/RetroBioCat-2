from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List


@dataclass
class TemplateFreeResult:
    """A single prediction from a template-free model."""
    reactants: str  # SMILES string of reactants (dot-separated)
    score: float  # model confidence / probability


class TemplateFreeModel(ABC):
    """Abstract base class for template-free retrosynthesis models.

    These models predict reactant SMILES directly from a product SMILES
    (no SMARTS templates). They use beam search to return multiple predictions.
    """

    @abstractmethod
    def predict(self, smi: str, num_results: int = 10) -> List[TemplateFreeResult]:
        """Predict reactants for a given product SMILES.

        Args:
            smi: Product SMILES string.
            num_results: Maximum number of predictions to return.

        Returns:
            List of TemplateFreeResult sorted by descending score.
        """
        pass

    @abstractmethod
    def is_loaded(self) -> bool:
        """Whether the model weights have been loaded."""
        pass

    @abstractmethod
    def load(self):
        """Load model weights. Called lazily on first predict() call."""
        pass
