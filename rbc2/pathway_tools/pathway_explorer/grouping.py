from typing import List

from rbc2 import Pathway

SCORES_TO_USE = ['change_in_complexity',
                 'fraction_starting_material_available',
                 'num_enzyme_steps_with_precedent']


def _get_scores_list(pathway: Pathway):
    scores_list = []
    scores_dict = pathway.scores
    for name in scores_dict:
        if name in SCORES_TO_USE:
            scores_list.append(scores_dict[name])
    return scores_list


def _generate_end_nodes_dict(pathways: List[Pathway]):
    """
    Iterates over pathways, sorting pathways into groups of: [end_nodes][reactions][scores]

    Pathways which are identical for these three are grouped together.
    """

    end_nodes_dict = {}
    for pathway in pathways:
        end_nodes = str(sorted(pathway.end_smis()))
        if end_nodes not in end_nodes_dict:
            end_nodes_dict[end_nodes] = {}

        reactions = [list(reaction.template_metadata.keys())[0] for reaction in pathway.reactions]
        reactions = str(sorted(reactions))
        if reactions not in end_nodes_dict[end_nodes]:
            end_nodes_dict[end_nodes][reactions] = {}

        scores = _get_scores_list(pathway)
        scores = str(scores)
        if scores not in end_nodes_dict[end_nodes][reactions]:
            end_nodes_dict[end_nodes][reactions][scores] = []
        end_nodes_dict[end_nodes][reactions][scores].append(pathway)

    return end_nodes_dict


def _get_grouped_pathways(end_nodes_dict):
    """ From end_nodes_dict, makes a list of grouped_pathway lists"""
    grouped_pathways = []
    for end_nodes in end_nodes_dict:
        for reactions in end_nodes_dict[end_nodes]:
            for scores in end_nodes_dict[end_nodes][reactions]:
                grouped_pathways.append(end_nodes_dict[end_nodes][reactions][scores])

    return grouped_pathways


def _collapse_groups(groups):
    """
    From a group of pathways, choose one to represent the rest
    """

    def sort_group(group):
        scores = []
        pathway_indexes = []
        for i, pathway in enumerate(group):
            enzyme_score = float(pathway.scores['num_enzyme_steps_with_precedent'])
            scores.append(enzyme_score)
            pathway_indexes.append(i)

        zipped_pairs = zip(scores, pathway_indexes)

        sorted_indexes = [x for _, x in sorted(zipped_pairs)]

        sorted_group = []
        for index in sorted_indexes:
            sorted_group.append(group[index])

        return sorted_group

    pathways = []

    for group in groups:
        sorted_group = sort_group(group)
        pathway = sorted_group.pop(0)
        for other_pathway in sorted_group:
            pathway.other_varients.append(other_pathway)
            pathway.other_varients_as_nodes.append(other_pathway.list_nodes)
        pathways.append(pathway)

    return pathways


def group_pathways(pathways: List[Pathway]) -> List[List[Pathway]]:
    end_nodes_dict = _generate_end_nodes_dict(pathways)  # groups by end_nodes, reactions, scores
    grouped_pathways = _get_grouped_pathways(end_nodes_dict)  # converts dict to a list of lists
    #new_pathways = _collapse_groups(grouped_pathways)
    return grouped_pathways
