import time
from typing import Optional, List

from rbc2.expansion.multi_expander import MultiExpander
from rbc2.reaction_evaluation.feasability import Filter, default_filter_repo
from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator import \
    DefaultSQLStartingMaterialEvaluator
from rbc2.utils.add_logger import add_logger
from rbc2.configs.logging_config import logging_config
from rbc2.configs.mcts_config import MCTS_Config
from rbc2.expansion.expander_repository import get_expanders
from rbc2.expansion.expander_interface import Expander
from rbc2.mcts.mcts_loop.backpropogate import backpropogate
from rbc2.mcts.mcts_loop.expansion.expand import Expansion
from rbc2.mcts.mcts_loop.rollout import rollout
from rbc2.mcts.mcts_loop.score_node import score_node
from rbc2.mcts.mcts_loop.selection import Selection
from rbc2.mcts.tree_node import create_root, MCTS_Node
from rbc2.data_model.network import Network
from rbc2.data_model.pathway import Pathway

class MCTS():

    def __init__(self,
                 target_smi: str,
                 expanders: dict[str: Expander],
                 filters: dict[str: Filter] = default_filter_repo,
                 starting_material_evaluator: Optional[DefaultSQLStartingMaterialEvaluator] = None,
                 network: Optional[Network] = None,
                 mcts_config: Optional[MCTS_Config] = None):

        self.target_smi = target_smi
        self.logger = add_logger('MCTS', level=logging_config.mcts)

        # config
        self.mcts_config = mcts_config
        if self.mcts_config is None:
            self.mcts_config = MCTS_Config()

        # starting material evaluator
        self.starting_material_evaluator = starting_material_evaluator
        if self.starting_material_evaluator is None:
            self.starting_material_evaluator = DefaultSQLStartingMaterialEvaluator()

        # network - used to save expansions so they are only done once
        self.network = network
        if self.network is None:
            self.network = Network()

        # multi_expander made up of the individual expanders
        self.multi_expander = MultiExpander(expanders, network=self.network)

        # filters
        self.filters = filters

        # mcts steps
        self.selection = Selection()
        self.expansion = Expansion(self.multi_expander,
                                   self.starting_material_evaluator,
                                   self.mcts_config)


        self.root: MCTS_Node = create_root(target_smi)  # root node
        self.solved = []  # the solved nodes, updated during the backpropagation step

        self.search_complete = False  # used to stop the search either on max iterations or max run time

        # stats
        self.iterations = 0
        self.run_time = 0
        self.positive_backpropagations = 0

    def run(self, callback=None):
        """Runs the MCTS search"""
        self.logger.debug(f'Running MCTS search for {self.target_smi}.  Max time: {self.mcts_config.max_search_time} seconds.  Max iterations: {self.mcts_config.max_iterations}')
        t0 = time.time()
        while self.search_complete is False:
            self.do_a_loop()
            self._check_run_time(t0)
            if callback is not None and self.iterations % self.mcts_config.callback_iterations == 0:
                callback(self)

    def do_a_loop(self):
        self.logger.debug(f'---- ITERATION {self.iterations} ----')
        node = self.selection.select(self.root, self.mcts_config.exploration)
        new_node = rollout(node, self.expansion, self.selection, self.network, self.filters, self.mcts_config)
        if new_node is None:
            self.logger.debug(f'Search complete - fully explored')
            self.search_complete = True
        score = score_node(new_node, self.mcts_config, self.starting_material_evaluator)
        if score >= 0.95: self.positive_backpropagations += 1
        self.solved += backpropogate(new_node, score)
        self.iterations += 1

    def _get_nodes(self, node: MCTS_Node) -> List[MCTS_Node]:
        """Returns all the nodes which are decendents of the given node"""
        nodes = []
        evaluated_children = [child for child in node.children if child.is_evaluated()]
        nodes += evaluated_children
        for child in evaluated_children:
            nodes += self._get_nodes(child)
        return nodes

    def get_all_nodes(self) -> List[MCTS_Node]:
        nodes = [self.root]
        nodes += self._get_nodes(self.root)
        return nodes

    def get_solved_nodes(self) -> List[MCTS_Node]:
        return self.solved

    def get_solved_pathways(self) -> List[Pathway]:
        solved = self.get_solved_nodes()
        pathways = [node.pathway for node in solved]
        return list(set(pathways))

    def _check_run_time(self, t0):
        if self.mcts_config.max_iterations is not None:
            if self.iterations >= self.mcts_config.max_iterations:
                self.logger.debug(f'Search complete after max iterations {self.iterations}')
                self.search_complete = True
                return

        self.run_time = round(time.time() - t0, 2)
        if self.mcts_config.max_search_time is not None:
            if self.run_time > self.mcts_config.max_search_time:
                self.logger.debug(f'Search complete after max time {round(self.run_time, 2)} seconds')
                self.search_complete = True
                return

    def get_run_stats(self):
        return {
            'target_smi': self.target_smi,
            'selections': self.selection.metrics,
            'template_applications': self.multi_expander.template_application_counts(),
            'expander_calls': self.multi_expander.expander_calls(),
            'iterations': self.iterations,
            'positive_backpropagations': self.positive_backpropagations,
            'run_time': self.run_time,
            'reactions': len(self.network.all_reactions()),
            'molecules': len(self.network.all_smis()),
            'solved': len(self.solved)
        }

    def get_flat_run_stats(self):
        stats = self.get_run_stats()
        flat_stats = {'target_smi': stats['target_smi'],
                      'iterations': stats['iterations'],
                      'positive_backpropagations': stats['positive_backpropagations'],
                      'run_time': stats['run_time'],
                      'reactions': stats['reactions'],
                      'molecules': stats['molecules'],
                      'solved': stats['solved']}
        for key, value in stats['selections'].items():
            flat_stats[f"selections_{key}"] = value
        for key, value in stats['template_applications'].items():
            flat_stats[f"template_applications_{key}"] = value
        for key, value in stats['expander_calls'].items():
            flat_stats[f"expander_calls_{key}"] = value
        return flat_stats

    def get_pa_routes(self):
        """ Outputs pathways in the pa routes style"""
        return [pathway.get_pa_route(self.starting_material_evaluator) for pathway in self.get_solved_pathways()]


if __name__ == '__main__':
    target_smi = 'C[C@H](C[C@H](C[C@@H](C[C@H](C[C@H](C)C(=O)O)C(=O)O)C[C@H](O)[C@H](O)C[C@@H](C[C@H](C[C@H](C)C(=O)O)C(=O)O)C[C@@H](C[C@@H](C)C(=O)O)C(=O)O)C(=O)O)C(=O)O'
    expanders = get_expanders(['retrobiocat', 'aizynthfinder'])
    mcts = MCTS(target_smi, expanders)
    mcts.mcts_config.max_length = 5
    mcts.mcts_config.max_search_time = 60*20  # max 20 minutes
    mcts.mcts_config.max_iterations = 1250
    mcts.mcts_config.chemistry_filter = 'None'
    mcts.multi_expander.expander_config.precedent_search = False
    mcts.run()

