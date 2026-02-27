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
                     rxn_domain: str = '',
                     fwd_rxn=False,
                     atom_mappings: Optional[dict] = None) -> List[Reaction]:

    if metadata is None:
        metadata = {}

    new_reactions = []
    for rxn_name, outcome_list in outcomes_dict.items():
        rxn_metadata = metadata.get(rxn_name, {})
        mapping_list = atom_mappings.get(rxn_name, []) if atom_mappings else []

        for i, rxn_outcome in enumerate(outcome_list):
            unique_id = make_unique_id(rxn_name)

            if fwd_rxn == True:
                # then the target smi is the substrate(s), and the outcomes are products.  Flip these for creation.
                product = '.'.join(rxn_outcome)
                substrates = target_smi.split('.')  # list of substrates is expected
            else:
                product = target_smi
                substrates = rxn_outcome

            data = {}
            if atom_mappings is not None:
                mapping = mapping_list[i] if i < len(mapping_list) else {}
                if isinstance(mapping, list):
                    # Multi-step: already a chain of step dicts
                    data['atom_mapping'] = mapping
                else:
                    # Single-step: wrap flat mapping dict in step format
                    data['atom_mapping'] = [
                        {
                            'rxn_name': rxn_name,
                            'target_smi': product,
                            'substrate_smis': substrates,
                            'mapping': mapping
                        }
                    ]

            reaction = Reaction(product=product,
                                substrates=substrates,
                                unique_id=unique_id,
                                name=rxn_name,
                                rxn_type=rxn_type,
                                rxn_domain=rxn_domain,
                                template_metadata={rxn_name: rxn_metadata},
                                score=score,
                                data=data)
            new_reactions.append(reaction)

    return new_reactions




