from rbc2.mcts.mcts import MCTS
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout

def to_networkx(mcts: MCTS):

    graph = nx.DiGraph()
    nodes = mcts.get_all_nodes()
    graph.add_node(mcts.root.id)
    for node in nodes:
        evaluated_children = [child for child in node.children if child.is_evaluated()]
        for child in evaluated_children:
            graph.add_edge(node.id, child.id)

    for node in nodes:
        if node.terminal == True:
            graph.nodes[node.id]['colour'] = 'red'
        if node.solved == True:
            graph.nodes[node.id]['colour'] = 'green'

    return graph


def draw_mcts(mcts):
    graph = to_networkx(mcts)
    pos = graphviz_layout(graph, prog="dot")

    graph_nodes = list(graph.nodes)
    graph_colours = []
    for node in graph_nodes:
        graph_colours.append(graph.nodes[node].get('colour', 'blue'))

    nx.draw(graph, pos, nodelist=graph_nodes, node_color=graph_colours, node_size=8, alpha=0.5, with_labels=False)
    plt.show()

if __name__ == '__main__':
    target_smi = '[C@H]1(C2=CC=CC=C2)NCCCC1'
    mcts = MCTS(target_smi, ('aizynthfinder',))
    mcts.config.max_length = 5
    mcts.config.chemistry_filter = None
    mcts.config.callback_iterations=50
    mcts.config.max_search_time = 120
    mcts.run(callback=draw_mcts)