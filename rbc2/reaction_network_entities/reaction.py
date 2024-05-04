import uuid
from dataclasses import dataclass, field, asdict
from typing import List, Optional
from rbc2.reaction_evaluation.complexity import get_complexity
from rbc2.reaction_network_entities.precedents import Precedent


@dataclass
class Reaction():
    product: str
    substrates: List[str]
    unique_id: str = ''
    name: str = 'unnamed_reaction'
    rxn_type: str = 'no_rxn_type'
    rxn_domain: str = 'no_rxn_domain'
    score: float = 0

    template_metadata: dict[str: dict] = field(default_factory=dict)
    precedents: List[Precedent] = field(default_factory=list)
    feasability_filter_scores: dict[str: float] = field(default_factory=dict)
    complexity_change: Optional[float] = None

    def __post_init__(self):
        if self.unique_id == '':
            self.unique_id = str(uuid.uuid4())

    def __hash__(self):
        return hash(self.unique_id)

    def __str__(self):
        return f'Reaction ({self.rxn_type}): {self.reaction_smiles()}'

    def get_complexity_change(self) -> float:
        if self.complexity_change is None:  # eg if product=4 and substrates=3, complexity_change=1
            self.complexity_change = round(
                get_complexity(self.product) - max([get_complexity(smi) for smi in self.substrates]), 3)
        return self.complexity_change

    def get_similarity_score(self) -> float:
        if len(self.precedents) == 0:
            return 0
        return self.precedents[0].similarity

    def reaction_smiles(self) -> str:
        rxn_smi = ".".join(self.substrates)
        rxn_smi += '>>'
        rxn_smi += self.product
        return rxn_smi

    def to_dict(self) -> dict:
        self.get_complexity_change()
        self.score = round(self.score, 3)
        return asdict(self)



def reaction_from_dict(reaction_dict) -> Reaction:
    """ Loads a reaction from a dictionary, such as would be returned by asdict(Reaction) """

    return Reaction(product=reaction_dict['product'],
                    substrates=reaction_dict['substrates'],
                    name=reaction_dict['name'],
                    rxn_type=reaction_dict['rxn_type'],
                    rxn_domain=reaction_dict['rxn_domain'],
                    unique_id=reaction_dict['unique_id'],

                    # these are optional, will just be empty/0 if not present
                    score=reaction_dict.get('score', 0),
                    precedents=[Precedent(**precedent_dict) for precedent_dict in reaction_dict.get('precedents', [])],
                    template_metadata=reaction_dict.get('template_metadata', {}),
                    feasability_filter_scores=reaction_dict.get('feasability_filter_scores', {}),
                    complexity_change=reaction_dict.get('complexity_change', None))

def reactions_to_dicts(reactions: List[Reaction]) -> List[dict]:
    """ Converts a list of reactions to a list of dictionaries, such as would be returned by asdict(Reaction) """

    rxn_dicts = [r.to_dict() for r in reactions]

    # drop None values
    rxn_dicts = [{k: v for k, v in rxn_dict.items() if v is not None} for rxn_dict in rxn_dicts]

    return rxn_dicts


def sort_reactions_by_score(reactions: List[Reaction]) -> List[Reaction]:
    return sorted(reactions, key=lambda x: x.score, reverse=True)