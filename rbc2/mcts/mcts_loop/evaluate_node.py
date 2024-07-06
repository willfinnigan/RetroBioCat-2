from rbc2.reaction_evaluation.feasability import Filter
from rbc2.template_application.create_reactions_from_output.find_duplicates import remove_duplicates
from rbc2.utils.add_logger import add_logger
from rbc2.mcts.mcts_logging_config import logging_config
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.mcts.tree_node import MCTS_Node, create_node_with_pathway
from rbc2.data_model.network import Network
from rbc2.data_model.reaction import Reaction
from rbc2.data_model.pathway import Pathway

node_evaluation_logger = add_logger('MCTS_Node_Evaluation', level=logging_config.mcts_node_evaluation)

def resolve_unevaluated_mcts_node(node: MCTS_Node, network: Network, tree_node_store: dict, filters: dict[str: Filter], mcts_config: MCTS_Config):
    """
    Run option to get evaluated MCTS Node(s)
    Remove unavaluated node from parents children
    """
    parent = node.parent
    parent.children.remove(node)

    new_mcts_nodes = []
    reactions = get_reactions_from_node_option(node, network, filters, mcts_config)
    for reaction in reactions:
        new_pathway = Pathway(parent.pathway.reactions + [reaction])
        new_node = create_node_with_pathway(parent, new_pathway, node.value)
        new_mcts_nodes.append(new_node)

    # use existing mcts nodes if they already exist
    for i, new_node in enumerate(new_mcts_nodes):
        new_node_hash = hash(new_node)
        if new_node_hash in tree_node_store:
            new_mcts_nodes[i] = tree_node_store[new_node_hash]
        else:
            tree_node_store[new_node_hash] = new_node

    if len(new_mcts_nodes) == 0:  # if evaluates to nothing, return parent
        if len(node.parent.children) == 0:  # if the parent has no children left, it is now terminal
            node.parent.terminal = True
        return node.parent
    else:
        node.parent.children += new_mcts_nodes
        return new_mcts_nodes[0]

def get_reactions_from_node_option(node: MCTS_Node, network: Network, filters: dict[str: Filter], mcts_config: MCTS_Config):
    if node.option.evaluate() == False:
        node.option.evaluate()

    to_keep = []
    for reaction in node.option.reactions:
        if are_reaction_substrates_already_in_the_parent_pathway(node.parent.pathway, reaction):
            continue
        if is_reaction_in_block_lists(reaction, mcts_config):
            continue
        if mcts_config.merge_reactions_from_same_domain:
            reaction = merge_reaction_if_match_already_existing_in_network_from_same_domain(reaction, network)
        to_keep.append(reaction)


    chem_fil = mcts_config.chemistry_filter
    chem_thres = mcts_config.chemistry_filter_cutoff
    biocat_fil = mcts_config.biocatalysis_filter
    feasible = [reaction for reaction in to_keep if does_reaction_pass_feasability_filters(reaction, filters, chem_fil, chem_thres, biocat_fil) is True]

    return feasible

def merge_reaction_if_match_already_existing_in_network_from_same_domain(reaction: Reaction, network: Network):

    if reaction in network.reactions:
        return reaction  # if reaction already in the network then no matching is needed

    reactions = remove_duplicates([reaction],
                                  network,
                                  require_same_domain=True,
                                  update_match_metadata=True,
                                  also_return_matches=True)
    if len(reactions) != 1:
        raise Exception(f'Error on merging reactions by domain during MCTS expansion - too many reactions returned ({len(reactions)})')

    return reactions[0]

def are_reaction_substrates_already_in_the_parent_pathway(pathway: Pathway, reaction: Reaction):
    return any(elem in pathway.all_smis for elem in reaction.substrates)

def is_reaction_in_block_lists(reaction: Reaction, mcts_config: MCTS_Config):
    return (mcts_config.avoid_blocked_reactions == True) and (reaction.unique_id in mcts_config.blocked_reactions)

def does_reaction_pass_feasability_filters(reaction: Reaction,
                                           filters: dict[str, Filter],
                                           chemistry_filter: str = 'None',
                                           chemistry_filter_cutoff: float = 0.05,
                                           biocatalysis_filter: str = 'None'):

    if reaction.rxn_domain == 'chemistry' and chemistry_filter != 'None':
        if chemistry_filter not in filters:
            node_evaluation_logger.warning(f'Chemistry filter {chemistry_filter} not found in filters dict')
            return True
        filter_func = filters[chemistry_filter]
        feasability = filter_func(reaction)
        return feasability >= chemistry_filter_cutoff

    if reaction.rxn_domain == 'biocatalysis' and biocatalysis_filter != 'None':
        if biocatalysis_filter not in filters:
            node_evaluation_logger.warning(f'Biocatalysis filter {biocatalysis_filter} not found in filters dict')
            return True
        filter_func = filters[biocatalysis_filter]
        feasability = filter_func(reaction)
        return feasability > 0

    return True