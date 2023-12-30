class LoggingConfig():

    def __init__(self):
        self.mcts = 'INFO'
        self.mcts_log_every_n_loops = 20
        self.mcts_selection = 'WARNING'
        self.mcts_backpropagate = 'WARNING'
        self.mcts_rollout = 'WARNING'
        self.mcts_tree_node = 'WARNING'
        self.mcts_unevaluated_treenode = 'WARNING'
        self.mcts_node_evaluation = 'WARNING'
        self.mcts_feasability_filter = 'WARNING'
        self.mcts_scorer = 'WARNING'
        self.mcts_expander = 'WARNING'
        self.mcts_expander_criteria_checker = 'WARNING'
        self.mcts_node_evaluation = 'WARNING'


    def set_global_mode(self, global_mode):
        self.mcts = global_mode
        self.mcts_selection = global_mode
        self.mcts_expander = global_mode
        self.mcts_backpropagate = global_mode
        self.mcts_rollout = global_mode
        self.mcts_tree_node = global_mode
        self.mcts_unevaluated_treenode = global_mode
        self.mcts_node_evaluation = global_mode
        self.mcts_feasability_filter = global_mode
        self.mcts_scorer = global_mode

    def set_mcts_loops_only(self):
        self.mcts = 'DEBUG'
        self.mcts_tree_node = 'DEBUG'
        #self.mcts_rollout = 'DEBUG'
        #self.mcts_expander = 'DEBUG'
        #self.mcts_node_evaluation = 'DEBUG'
        #self.mcts_expander_criteria_checker = 'DEBUG'
        #self.mcts_selection = 'DEBUG'
        #self.mcts_backpropagate = 'DEBUG'
        #self.mcts_feasability_filter = 'DEBUG'


logging_config = LoggingConfig()
logging_config.set_global_mode('DEBUG')
#logging_config.mcts_expander = 'DEBUG'
#logging_config.mcts_selection = 'DEBUG'
#logging_config.mcts_rollout = 'DEBUG'
#logging_config.set_global_mode('DEBUG')
#logging_config.set_mcts_loops_only()