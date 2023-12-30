

class PathwayScoring_Config():

    def __init__(self):

        # clustering
        self.dont_cluster_single_step_routes = True
        self.distance_threshold = 2

        # costing
        self.in_stock_cost = 1
        self.not_in_stock_cost = 5
        self.reaction_cost = 0.5
        self.reaction_yield = 0.9

        # adjust non-buyable with complexity
        self.adjust_with_complexity = True
        self.complexity_multiplier = 2

    def update_from_dict(self, attr_dict):
        current_dict = self.to_dict()
        for key, value in attr_dict.items():
            if key in current_dict:
                setattr(self, key, value)
        return self

    def to_dict(self):
        return self.__dict__

def pathway_scoring_config_from_dict(config_dict: dict) -> PathwayScoring_Config:
    return PathwayScoring_Config().update_from_dict(config_dict)
