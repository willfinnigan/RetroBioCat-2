from collections.abc import Callable
from typing import List, Optional

from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluatorInterface
from rbc2.reaction_network_entities.reaction import Reaction
from rbc2.utils.rdkit_utils import rdkit_smile

""" Functions for converting pa_routes to and from other formats. """

def get_pa_route(smi: str,
                 starting_material_evaluator: StartingMaterialEvaluatorInterface,
                 get_smi_produced_by: Callable[[str], List[Reaction]]) -> dict:
    """
    Recursive function which will generate a 'pa_route' - format commonly used by AIZynthfinder associated tools.
    """

    in_stock = bool(starting_material_evaluator.eval(smi)[0])
    tree = {'type': 'mol', 'smiles': smi, 'children': [], 'in_stock': in_stock}
    for reaction in get_smi_produced_by(smi):
        rxn_branch = {'type': 'reaction', 'smiles': reaction.reaction_smiles(), 'children': []}
        for child_smi in reaction.substrates:
            childs_pa_route = get_pa_route(child_smi, starting_material_evaluator, get_smi_produced_by)
            rxn_branch['children'].append(childs_pa_route)

        if len(rxn_branch['children']) == 0:  # ensure no empty children lists
            rxn_branch.pop('children')
        tree['children'].append(rxn_branch)

    if len(tree['children']) == 0:  # ensure no empty children lists
        tree.pop('children')

    return tree

def get_route_target(pa_route: dict) -> str:
    """Returns the target smiles of a pa_route"""
    if pa_route['type'] == 'mol':
        return pa_route['smiles']
    else:
        raise Exception('Route target is not a molecule')


def get_route_leaves(pa_route: dict, leafs=None) -> List[str]:
    """Return a list of smiles for the ends of the pa_route"""
    if leafs is None: leafs = []

    if 'children' in pa_route:
        for child in pa_route['children']:
            leafs = get_route_leaves(child, leafs=leafs)
    elif (pa_route['type'] == 'mol'):
        leafs.append(pa_route['smiles'])

    return leafs

def get_all_smis(pa_route: dict, all_smis=None) -> List[str]:
    """Return a list of smiles for all the nodes in the pa_route"""
    if all_smis is None: all_smis = []

    if pa_route['type'] != 'mol':
        raise Exception('type is not mol')

    if pa_route.get('smiles', None) is not None:
        all_smis.append(pa_route['smiles'])
        for child_rxn in pa_route.get('children', []):
                for child_mol in child_rxn.get('children', []):
                    all_smis = get_all_smis(child_mol, all_smis=all_smis)

    return all_smis

def get_route_intermediates(pa_route: dict) -> List[str]:
    """Return a list of smiles which are not the target and not the leaves"""
    leaves = get_route_leaves(pa_route)
    target = get_route_target(pa_route)
    all_smis = get_all_smis(pa_route)
    intermediates = [smi for smi in all_smis if smi not in leaves and smi != target]
    return intermediates


def reactions_from_pa_route(pa_route: dict, reactions: Optional[List[Reaction]] = None, force_rdkit=True, rxn_name='test_set') -> List[Reaction]:
    """Return a list of reactions from a pa_route"""
    if reactions is None: reactions = []

    if pa_route['type'] == 'reaction':
        rxn_smi = pa_route['smiles']
        substrates = rxn_smi.split('>>')[0].split('.')
        product = rxn_smi.split('>>')[1]

        if force_rdkit == True:
            substrates = [rdkit_smile(smi) for smi in substrates]
            product = rdkit_smile(product)

        rxn = Reaction(substrates=substrates,
                       product=product,
                       name=rxn_name,
                       rxn_type=rxn_name,
                       rxn_domain=rxn_name)
        reactions.append(rxn)

    for child in pa_route.get('children', []):
        reactions += reactions_from_pa_route(child)

    return reactions


