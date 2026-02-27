from typing import Optional, List

from rdkit import Chem

from rbc2.configs.expansion_config import Expansion_Config
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction, sort_reactions_by_score
from rbc2.data_model.reaction_option import ReactionOption, sort_options_by_score
from rbc2.expansion.default_expander import DefaultExpander
from rbc2.expansion.expanders.template_free.models.model_interface import TemplateFreeModel, TemplateFreeResult
from rbc2.template_application.create_reactions_from_output.process_reactions import reaction_sanity_check, \
    process_reactions


class TemplateFreeExpander(DefaultExpander):
    """Base class for template-free retrosynthesis expanders.

    Template-free models predict reactant SMILES directly from product SMILES
    (no SMARTS templates). This class follows the same pre-evaluated option pattern
    as RetroBioCatExpander: overrides _create_new_options() to run the model,
    creates ReactionOption objects with evaluated=True and reactions already populated.
    """

    def __init__(self,
                 model: TemplateFreeModel,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None,
                 num_results: int = 10):
        super().__init__(network=network, config=config)
        self.model = model
        self.num_results = num_results
        # rxn_type, rxn_domain set by subclasses

    def create_option(self, smi: str, name: str, smarts: list,
                      template_metadata: dict, score: float) -> ReactionOption:
        option = ReactionOption(smi,
                                name,
                                [],  # no smarts for template-free
                                template_metadata,
                                self.rxn_type,
                                self.rxn_domain,
                                score,
                                None,  # no evaluation function - already evaluated
                                precedent_search_function=None)
        return option

    def _create_new_options(self, smi: str) -> List[ReactionOption]:
        self.calls += 1

        # Lazy load model
        if not self.model.is_loaded():
            self.model.load()

        # Run model prediction
        predictions = self.model.predict(smi, num_results=self.num_results)

        # Process predictions into reactions
        reactions = self._process_predictions(smi, predictions)

        if self.network is not None:
            [self.network.add_reaction(reaction) for reaction in reactions]

        reactions = sort_reactions_by_score(reactions)

        if self.config.max_reactions is not None:
            reactions = reactions[:self.config.max_reactions]

        # Create pre-evaluated options
        options = []
        for i, reaction in enumerate(reactions):
            option = self.create_option(
                smi=smi,
                name=f"{self.rxn_type}_{i}",
                smarts=[],
                template_metadata={},
                score=reaction.score)
            option.evaluated = True
            option.reactions = [reaction]
            options.append(option)

        sorted_options = sort_options_by_score(options)
        if self.network is not None:
            self.network.bulk_add_options(smi, self.rxn_type, sorted_options)

        return sorted_options

    def _process_predictions(self, product_smi: str,
                             predictions: List[TemplateFreeResult]) -> List[Reaction]:
        """Convert model predictions to validated Reaction objects."""
        reactions = []
        seen_substrate_sets = set()
        canon_product = self._canonicalize(product_smi)
        if canon_product is None:
            return []

        for pred in predictions:
            substrates = self._validate_and_canonicalize(pred.reactants)
            if substrates is None:
                continue

            # Deduplicate by substrate set
            substrate_key = tuple(sorted(substrates))
            if substrate_key in seen_substrate_sets:
                continue
            seen_substrate_sets.add(substrate_key)

            # Don't include reactions where product appears in substrates
            if canon_product in substrates:
                continue

            # Clamp score to [0, 1] regardless of model-specific normalization
            clamped_score = max(0.0, min(1.0, pred.score))

            reaction = Reaction(
                product=product_smi,
                substrates=list(substrates),
                name=f"{self.rxn_type}",
                rxn_type=self.rxn_type,
                rxn_domain=self.rxn_domain,
                score=clamped_score,
                template_metadata={})
            reactions.append(reaction)

        # Run sanity check
        reactions = reaction_sanity_check(reactions)

        # Process through network if available
        if self.network is not None:
            reactions = process_reactions(reactions, self.network, self.config)

        return reactions

    def _validate_and_canonicalize(self, reactants_smi: str) -> Optional[List[str]]:
        """Validate and canonicalize a dot-separated reactants SMILES string.
        Returns None if any component is invalid."""
        parts = reactants_smi.split('.')
        canonical = []
        for part in parts:
            c = self._canonicalize(part)
            if c is None:
                return None
            canonical.append(c)
        return canonical

    def _canonicalize(self, smi: str) -> Optional[str]:
        """Canonicalize a SMILES string. Returns None if invalid."""
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    def number_of_rule_applications(self) -> int:
        return 0  # No rule applications for template-free models
