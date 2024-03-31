from __future__ import annotations

import uuid
from copy import copy
from dataclasses import dataclass, field
from typing import Optional

from rbc2.data_model.reaction_option import ReactionOption
from rbc2.data_model.pathway import Pathway


@dataclass
class MCTS_Node():

    parent: Optional[MCTS_Node] = None
    pathway: Optional[Pathway] = None
    option: Optional[ReactionOption] = None

    terminal: bool = False
    children: list = field(default_factory=list)
    visits: int = 1
    value: float = 0
    solved: bool = False
    fully_searched: bool = False
    expanded: bool = False
    depth: int = 0
    is_root: bool = False

    def __post_init__(self):
        self.id = str(uuid.uuid4())

    def __hash__(self):
        return hash(self.id)

    def is_evaluated(self):
        if self.pathway is None and self.option is None:
            raise Exception("MCTS must either have a pathway (evaluated), or an option and a parent (non_evaluated)")
        if self.option is not None:
            if self.parent is None:
                raise Exception("If node is initialised with a ReactionOption, it must have a parent node")
            elif self.parent.pathway is None:
                raise Exception("If node is initialised with a ReactionOption, it's parent node must have a pathway")

        return self.pathway is not None

    def get_last_rxn_type(self):
        if self.option is not None:
            return self.option.rxn_type
        if self.pathway is not None:
            if len(self.pathway.reactions) > 0:
                rxn_type = self.pathway.reactions[-1].rxn_type
                return rxn_type
            elif self.is_root == False:
                raise Exception("Pathway has no reactions, but is not root")
        return None

def create_root(target_smi: str) -> MCTS_Node:
    pathway = Pathway([], target_smi=target_smi)
    root = MCTS_Node(pathway=pathway, is_root=True)
    return root

def create_node_from_option(parent: MCTS_Node, option: ReactionOption) -> MCTS_Node:
    depth = parent.depth + 1
    return MCTS_Node(parent=parent, option=option, value=copy(option.score), depth=depth)

def create_node_with_pathway(parent: MCTS_Node, pathway: Pathway, value: float) -> MCTS_Node:
    depth = parent.depth + 1
    return MCTS_Node(parent=parent, pathway=pathway, value=value, depth=depth)