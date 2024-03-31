import uuid
from typing import List, Optional

from rbc2.data_model.reaction import Reaction

# outcomes_dict will = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1'], ['C1=N[C@H](c2ccccc2)CCC1']]}

ID_SPLIT_CHAR = "_id_"

def make_unique_id(rxn_name):
    return f"{rxn_name}{ID_SPLIT_CHAR}{str(uuid.uuid4())}"

def create_reactions(target_smi: str,
                     outcomes_dict: dict,
                     score: float = 0,
                     metadata: Optional[dict] = None,
                     rxn_type: str = '',
                     rxn_domain: str = '') -> List[Reaction]:

    if metadata is None:
        metadata = {}

    new_reactions = []
    for rxn_name, outcome_list in outcomes_dict.items():
        rxn_metadata = metadata.get(rxn_name, {})
        for rxn_outcome in outcome_list:
            unique_id = make_unique_id(rxn_name)
            reaction = Reaction(product=target_smi,
                                substrates=rxn_outcome,
                                unique_id=unique_id,
                                name=rxn_name,
                                rxn_type=rxn_type,
                                rxn_domain=rxn_domain,
                                template_metadata={rxn_name: rxn_metadata},
                                score=score)
            new_reactions.append(reaction)

    return new_reactions

