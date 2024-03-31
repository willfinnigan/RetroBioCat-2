from typing import List
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator import \
    DefaultSQLStartingMaterialEvaluator
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluatorInterface
from rbc2.data_model.pathway import Pathway
from rbc2.data_model.reaction import Reaction


def leaf_molecule_availability(pathway: Pathway,
                               starting_material_evaluator: StartingMaterialEvaluatorInterface):

    available, not_available = [], []
    target_smi = pathway.target_smi
    list_end_smis = list(pathway.end_smi_depths.keys())

    if target_smi in list_end_smis and starting_material_evaluator.config.target_always_not_buyable is True:
        list_end_smis.remove(target_smi)
        not_available.append(target_smi)

    for smi in list_end_smis:
        in_stock, in_stock_info = starting_material_evaluator.eval(smi)
        if in_stock == True:
            available.append(smi)
        else:
            not_available.append(smi)

    return available, not_available


IN_STOCK_COST = 1
NOT_IN_STOCK_COST = 10
REACTION_YIELD = 0.8
REACTION_COST = 1

def cost_pathway(pathway: Pathway, starting_material_evaluator: StartingMaterialEvaluatorInterface) -> float:
    node_costs = _get_end_node_costs(pathway, starting_material_evaluator)   # firstly cost the end nodes by in stock or not
    nodes = pathway.end_smis()  # initially nodes is just the leaves
    while pathway.target_smi not in nodes:
        reactions = _get_reactions_of_smis(nodes, pathway)
        nodes = []
        for reaction in reactions:
            if _are_all_smis_costed(reaction.substrates, node_costs):
                node_costs[reaction.product] = _cost_intermediate(reaction.substrates, node_costs)
                nodes.append(reaction.product)

    score = float(node_costs[pathway.target_smi])
    return score


def rank_pathways(pathways, scores):
    """Ranks pathways by their score

    Args:
        pathways (list): list of pathways
        scores (list): list of scores

    Returns:
        list: ranked pathways
    """
    return [p for _, p in sorted(zip(scores, pathways), reverse=False)]


def _cost_intermediate(substrate_smis: List[str], node_costs: dict) -> float:
    cost_precursors = 0
    for smi in substrate_smis:
        cost_precursors += node_costs[smi]

    cost = (cost_precursors / REACTION_YIELD) + REACTION_COST
    return cost


def _get_reactions_of_smis(list_smis: List[str], pathway) -> List[Reaction]:
    reactions = []
    for smi in list_smis:
        reactions += list(pathway.smi_substrate_of[smi])
    return reactions


def _are_all_smis_costed(list_smis: List[str], node_costs: dict[str: float]) -> bool:
    for smi in list_smis:
        if node_costs.get(smi, None) is None:
            return False
    return True


def _get_end_node_costs(pathway: Pathway, starting_material_evaluator: StartingMaterialEvaluatorInterface) -> dict[str: float]:
    end_node_costs = {}
    buyable, not_buyable = leaf_molecule_availability(pathway, starting_material_evaluator)
    for smi in buyable:
        end_node_costs[smi] = IN_STOCK_COST
    for smi in not_buyable:
        end_node_costs[smi] = NOT_IN_STOCK_COST
        #if config.adjust_with_complexity:
            #mol = pathway.mol_dict[smi]
            #end_node_costs[smi] += mol.get_complexity() * float(self.config.complexity_multiplier)
    return end_node_costs