from collections import defaultdict
from dataclasses import asdict
from typing import List, Optional

from rbc2.pathway_tools.pa_route_conversion import get_pa_route
from rbc2.reaction_evaluation.starting_material_evaluator import StartingMaterialEvaluator
from rbc2.reaction_network_entities.reaction import Reaction, reaction_from_dict
from rbc2.utils.add_logger import add_logger

pathway_logger = add_logger('Pathway')

class Pathway:

    def __init__(self, reactions: List[Reaction], target_smi: Optional[str] = None):
        self.reactions = reactions

        self.smi_produced_by = defaultdict(set)
        self.smi_substrate_of = defaultdict(set)

        for reaction in self.reactions:
            self.smi_produced_by[reaction.product].add(reaction)
            for smi in reaction.substrates:
                self.smi_substrate_of[smi].add(reaction)

        self.product_smis = set(self.smi_produced_by.keys())
        self.substrate_smis = set(self.smi_substrate_of.keys())
        self.all_smis = self.product_smis | self.substrate_smis

        if target_smi is not None:
            self.target_smi = target_smi
            self.all_smis.add(self.target_smi)
        else:
            self.target_smi = self._get_target_smi()

        self.pathway_length = 0
        self.end_smi_depths: dict[str: int] = {}
        self.tree = self._make_tree(self.target_smi)

    def _get_target_smi(self):
        target_smis = [smi for smi in self.product_smis if smi not in self.substrate_smis]
        if len(target_smis) > 1:
            raise Exception('Pathway has multiple targets')
        elif len(target_smis) == 0:
            raise Exception('Pathway has no target')
        return target_smis[0]

    def _make_tree(self, smi: str, depth=0) -> dict:
        if self.pathway_length < depth:
            self.pathway_length = depth

        tree = {'smiles': smi, 'depth': depth, 'children': []}
        for reaction in self.smi_produced_by[smi]:
            for child_smi in reaction.substrates:
                tree['children'].append(self._make_tree(child_smi, depth=depth+1))

        if len(self.smi_produced_by[smi]) == 0:
            self.end_smi_depths[smi] = depth

        return tree

    def get_pa_route(self, starting_material_evaluator: StartingMaterialEvaluator):
        def get_smi_produced_by(smi):
            return list(self.smi_produced_by[smi])

        return get_pa_route(self.target_smi, starting_material_evaluator, get_smi_produced_by)

    def get_smi_producted_by(self, smi: str) -> Reaction:
        reactions = self.smi_produced_by[smi]
        if len(reactions) != 1:
            raise Exception(f'smi {smi} produced by multiple reactions')
        return list(reactions)[0]

    def end_smis(self):
        return list(self.end_smi_depths.keys())

    def save(self):
        """Returns a list of dicts containing the reactions in the pathway"""
        return [asdict(reaction) for reaction in self.reactions]

    def get_reaction_with_product(self, smi: str) -> Optional[Reaction]:
        reactions = self.smi_produced_by[smi]
        if len(reactions) == 0:
            return None

        if len(reactions) != 1:
            pathway_logger.warning(f'smi {smi} produced by multiple reactions')

        return list(reactions)[0]

    def get_reaction_with_substrate(self, smi: str) -> Optional[Reaction]:
        reactions = self.smi_substrate_of[smi]

        if len(reactions) == 0:
            return None

        if len(reactions) != 1:
            pathway_logger.warning(f'smi {smi} substrate of multiple reactions')

        return list(reactions)[0]



def load_pathway(reaction_dict_list: List[dict]):
    """Loads a list of dicts containing the reactions in the pathway"""
    return Pathway([reaction_from_dict(reaction_dict) for reaction_dict in reaction_dict_list])