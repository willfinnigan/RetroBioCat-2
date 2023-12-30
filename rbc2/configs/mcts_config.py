from rbc2.reaction_evaluation.feasability_filters import AIZYNTHFINDER_FILTER


class MCTS_Config():

    def __init__(self):

        # mcts setting
        self.max_length = 4
        self.exploration = 1.414  # optimal exploration value for UCB1 is sqrt(2)=1.414 if score [0-1]
        self.max_search_time = 120
        self.max_iterations = None
        self.callback_iterations = 20  # number of iterations before the mcts callback function is called (if set)

        # mcts scoring
        self.use_reaction_scores_for_mcts_initial_values = True  # recommended this is True, otherwise first round will be random
        self.score_mode = 'basic'  # 'basic', 'complexity_penalty', 'mass_percent', ('number_of_atoms'-not implemented yet)
        self.use_pathway_length_score = True  # also use a pathway length score (like aizynthfinder)

        # values for complexity penalty (if used)
        self.non_buyable_score = 0.2  # the default score for a non buyable compound
        self.max_complexity_penalty = -0.2  # the maximum penalty for a *complex* non buyable compound
        self.rel_complexity_no_penalty = 0  # complexity above this has no penalty
        self.rel_complexity_max_penalty = -1  # complexity below this has max penalty

        # multi_expansion options
        self.option_combination_method = 'order_by_score'  # ['interleave, order_by_score']

        # expansion option scoring
        self.allow_moves_beyond_solved = 0  # number of moves beyond a solved node that are allowed, normally 0
        self.stop_expansion_if_nonbuyable_at_max_length = False  # dont expand a mcts_node if a non buyable is found at max length (its impossible to solve)
        self.boost_enzyme_score_if_in_cascade = False
        self.boost_enzyme_in_cascade_score_by = 0.2



        self.max_chemistry_nodes = None

        self.chemistry_only_at_beginning_or_end = False  # if true, only allow chemistry at beginning or end of pathway
        self.max_chemistry_at_beginning = None  # if only allowing chemistry at beginning or end, optionally set a max number of chemistry nodes at beginning
        self.max_chemistry_at_end = None # if only allowing chemistry at beginning or end, optionally set a max number of chemistry nodes at end

        # expansion node evaluation
        self.avoid_blocked_reactions = True
        self.blocked_reactions = []

        self.merge_reactions_from_same_domain = False

        self.chemistry_filter = AIZYNTHFINDER_FILTER
        self.chemistry_filter_cutoff = 0.05
        self.biocatalysis_filter = 'None'  # RETROBIOCAT_FILTER


    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

def mcts_config_from_dict(config_dict: dict) -> MCTS_Config:
    return MCTS_Config().update_from_dict(config_dict)
