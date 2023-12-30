from dataclasses import dataclass, field
from typing import List

import yaml

from rbc2.configs.data_path import path_to_data_folder
from rbc2.expansion.expanders.retrobiocat_reaction_retrieval.rxn_class_interface import RetroBioCatReactions, \
    rxns_dict, multi_rxns_dict
from rbc2.template_application.apply_template.rdchiral.initialization import rdchiralReaction
from rbc2.utils.add_logger import add_logger

data_folder = f'{path_to_data_folder}/retrobiocat'


@dataclass
class RetroBioCatReaction:
    name: str = 'no_name'
    smarts: List[str] = field(default_factory=list)
    enzymes: List[str] = field(default_factory=list)
    cofactors: dict = field(default_factory=dict)
    positive_tests: List[str] = field(default_factory=list)
    negative_tests: List[str] = field(default_factory=list)
    example_rxn_string: str = ''
    type: str = ''
    experimental: bool = False
    two_step: bool = False
    requires_absence_of_water: bool = False

    product_seeds: List[str] = field(default_factory=list)
    substrate_1_seeds: List[str] = field(default_factory=list)
    substrate_2_seeds: List[str] = field(default_factory=list)

    # reaction can have no smarts, but reference multiple other reactions which can be run in sequence
    steps: List[List[str]] = field(default_factory=list)


def load_reactions_from_yaml() -> List[RetroBioCatReaction]:
    with open(f"{data_folder}/rxns_yaml.yaml") as stream:
        data_loaded = yaml.safe_load(stream)

    reactions = []
    for reaction_name, values in data_loaded.items():
        values['name'] = reaction_name
        rxn = RetroBioCatReaction(**values)
        reactions.append(rxn)

    return reactions


class YAML_RetroBioCatReactions(RetroBioCatReactions):

    def __init__(self,
               include_experimental=False,
               include_two_step=True,
               include_requires_absence_of_water=False):

        self.include_experimental = include_experimental
        self.include_two_step = include_two_step
        self.include_requires_absence_of_water = include_requires_absence_of_water

        self.logger = add_logger('YAMLRetroBioCat', level='WARNING')
        self.rxns = {}
        self.multi_rxns = {}
        self.reactions, self.enzymes, self.reaction_enzyme_map = [], [], {}
        self.reactionEnzymeCofactorDict = {}
        self.rxns_strings = {}
        self.rules_by_type = {}

    def get_rxns(self) -> rxns_dict:
        if self.rxns == {}:
            self._load_rxns()
        return self.rxns

    def get_multistep_rxns(self) -> multi_rxns_dict:
        if self.multi_rxns == {}:
            self._load_rxns()
        return self.multi_rxns

    def get_possible_enzymes(self, rxn_name) -> list:
        return self.reaction_enzyme_map.get(rxn_name, [''])

    def get_rxn_cofactors(self, rxn_name) -> dict:
        return self.reactionEnzymeCofactorDict.get(rxn_name, {})

    def get_enzyme_full_name(self, enzyme) -> str:
        return enzyme


    def _load_rxns(self):
        if self.rxns != {}:
            return

        reactions = load_reactions_from_yaml()

        # filter by options
        if self.include_experimental == False:
            reactions = [r for r in reactions if r.experimental != True]
        if self.include_requires_absence_of_water == False:
            reactions = [r for r in reactions if r.requires_absence_of_water != True]
        if self.include_two_step == False:
            reactions = [r for r in reactions if r.two_step != True]

        reaction_names = set()
        enzyme_names = set()

        # load single step reactions
        for reaction in reactions:
            if reaction.two_step == True:
                continue

            self.rxns[reaction.name] = [rdchiralReaction(smarts) for smarts in reaction.smarts]
            reaction_names.add(reaction.name)

            for enzyme in reaction.enzymes:
                enzyme_names.add(enzyme)
                self.reaction_enzyme_map.setdefault(reaction.name, set()).add(enzyme)

            self.reactions = list(reaction_names)
            self.enzymes = list(enzyme_names)
            self.reaction_enzyme_map = {k: list(v) for k, v in self.reaction_enzyme_map.items()}

        # load multi_step reactions
        for reaction in reactions:
            if reaction.two_step == True and len(reaction.steps) != 0:
                # this will become a list of lists of lists: [[[r1a, r1b], [r2a, r2b]], [[r1a, r1b], [r2a, r2b]]]
                self.multi_rxns[reaction.name] = self._get_multi_smarts_from_named_steps(self.rxns, reaction.steps)


    @staticmethod
    def _get_multi_smarts_from_named_steps(rxns, steps):
        processed_steps = []
        for group_steps in steps:
            group_smas = []
            for step_rxn in group_steps:
                step = []
                step_rdchiral_rxns = rxns[step_rxn]
                for rdchiral_rxn in step_rdchiral_rxns:
                    step.append(rdchiral_rxn)
                group_smas.append(step)

            processed_steps.append(group_smas)

        return processed_steps


if __name__ == '__main__':

    rxn_class = YAML_RetroBioCatReactions()
    smi = 'CC(C)(C)OC(=O)N1CCC(CC1)N(C)C(=O)C1=CC=C(C=C1)C1=CC=CC=C1'
    rxns = rxn_class.get_rxns()
    multi_step_rxns = rxn_class.get_multistep_rxns()

    for rxn, smarts in rxns.items():
        print(f"{rxn} - {smarts}")

    for multi_rxn, multi_smarts in multi_step_rxns.items():
        print(f"{multi_rxn} - {multi_smarts}")