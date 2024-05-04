from __future__ import annotations

from collections import defaultdict
from typing import List, Set,  Sequence

from typing import TYPE_CHECKING

from rbc2.pathway_tools.pa_route_conversion import get_pa_route

from rbc2.reaction_evaluation.starting_material_evaluator.starting_material_evaluator_interface import \
    StartingMaterialEvaluatorInterface

from rbc2.data_model.reaction import reactions_to_dicts, reaction_from_dict

if TYPE_CHECKING:
    from rbc2.data_model.reaction import Reaction
    from rbc2.data_model.reaction_option import ReactionOption
    from rbc2.expansion.expander_interface import Expander

ReactionID = str
OptionID = str
ExpanderID = str
Smi = str
RxnType = str

class Network():
    """ Network is used to keep a record of the outcome of all expansions."""

    def __init__(self, reactions: Sequence[Reaction] = ()):

        self.smi_produced_by: dict[Smi: Set[Reaction]] = defaultdict(set)
        self.smi_substrate_of: dict[Smi: Set[Reaction]] = defaultdict(set)
        self.reaction_options: dict[Smi: dict[ExpanderID: List[ReactionOption]]] = defaultdict(lambda: defaultdict(dict))
        self.reactions: Set[Reaction] = set()

        if len(reactions) != 0:
            for rxn in reactions:
                self.add_reaction(rxn)

    def add_reaction(self, reaction: Reaction):
        self.reactions.add(reaction)
        self.smi_produced_by[reaction.product].add(reaction)
        for smi in reaction.substrates:
            self.smi_substrate_of[smi].add(reaction)

    def remove_reaction(self, reaction: Reaction):
        self.reactions.discard(reaction)
        self.smi_produced_by[reaction.product].discard(reaction)
        for smi in reaction.substrates:
            self.smi_substrate_of[smi].discard(reaction)

    def add_option(self, option: ReactionOption):
        self.reaction_options[option.target_smi][option.rxn_type][option.unique_id] = option

    def bulk_add_options(self, smi: Smi, rxn_type: RxnType, list_options: List[ReactionOption]):
        self.reaction_options[smi][rxn_type] = {option.unique_id: option for option in list_options}

    def remove_option(self, option: ReactionOption):
        self.reaction_options[option.target_smi][option.rxn_type].pop(option.unique_id, None)

    def get_reaction_options(self, smi: Smi, rxn_type: RxnType) -> list[ReactionOption]:
        options_for_smi = self.reaction_options.get(smi, {})
        options_for_rxn_type = options_for_smi.get(rxn_type, {})
        return list(options_for_rxn_type.values())

    def are_options_available(self, smi: Smi, rxn_type: RxnType) -> bool:
        return self.reaction_options.get(smi, {}).get(rxn_type, False) is not False

    def get_reactions_which_molecule_is_produced_by(self, smi: Smi) -> Set[Reaction]:
        return self.smi_produced_by.get(smi, set())

    def get_reactions_which_molecule_is_substrate_of(self, smi: Smi) -> Set[Reaction]:
        return self.smi_substrate_of.get(smi, set())

    def all_smis(self) -> Set[Smi]:
        all_smis = set(self.smi_produced_by.keys())
        all_smis.update(set(self.smi_substrate_of.keys()))
        return all_smis

    def all_reactions(self) -> List[Reaction]:
        return list(self.reactions)

    def all_reaction_options(self) -> List[ReactionOption]:
        all_options = []
        for smi, rxn_type_options in self.reaction_options.items():
            for rxn_type, options_dict in rxn_type_options.items():
                for option_id, option in options_dict.items():
                    all_options.append(option)
        return all_options

    def save(self):
        """Save the network to a dict"""
        data = {"reactions": reactions_to_dicts(self.all_reactions()),
                "reaction_options": [option_to_dict(opt) for opt in self.all_reaction_options()]}
        return data

    def load(self, data: dict, expanders: List[Expander]):
        """
        Load the network from data dict
        ReactionOptions will only be loaded if the relevant expander is provided
        """

        # check each expander is associated with this network
        for expander in expanders:
            if expander.network != self:
                raise Exception("Can not load reaction options when expander is not associated with the same network")

        # load reactions
        reaction_unique_id_dict = {}
        for reaction_dict in data['reactions']:
            reaction = reaction_from_dict(reaction_dict)
            reaction_unique_id_dict[reaction.unique_id] = reaction
            self.add_reaction(reaction)

        # load reaction options
        expander_dict = {exp.rxn_type: exp for exp in expanders}
        for opt_dict in data['reaction_options']:
            rxn_type = opt_dict['rxn_type']
            expander = expander_dict.get(rxn_type, None)
            if expander is None:
                continue

            option = option_from_dict(opt_dict, expander)

            # add reactions from ids
            for unique_id in opt_dict.get('reaction_ids', []):
                reaction = reaction_unique_id_dict.get(unique_id, None)
                if reaction is None:
                    continue
                option.reactions.append(reaction)

            self.add_option(option)


    def get_pa_route(self, start_smi, starting_material_evaluator: StartingMaterialEvaluatorInterface):
        def get_smi_produced_by(smi):
            return list(self.smi_produced_by[smi])
        return get_pa_route(start_smi, starting_material_evaluator, get_smi_produced_by)


def option_to_dict(option: ReactionOption) -> dict:
    """
    Saves the data of a ReactionOption to a dict
    Use the option_from_dict function to recreate the ReactionOption
    """
    opt_dict = {'target_smi': option.target_smi,
                'name': option.name,
                'smarts': option.smarts,
                'metadata': option.metadata,
                'rxn_type': option.rxn_type,
                'rxn_domain': option.rxn_domain,
                'score': option.score,
                'evaluated': option.evaluated,
                'precedents_searched': option.precedents_searched,
                'reaction_ids': [reaction.unique_id for reaction in option.reactions]
                }
    return opt_dict


def option_from_dict(option_dict: dict, expander: Expander) -> ReactionOption:
    """
    Load a previously saved option from a dict - using the expander method
    """

    if option_dict.get('rxn_type', '') != expander.rxn_type:
        raise Exception('Expander rxn_type does not match option_dict rxn_type')

    return expander.create_option(smi=option_dict['target_smi'],
                                  name=option_dict['name'],
                                  smarts=option_dict['smarts'],
                                  template_metadata=option_dict['metadata'],
                                  score=option_dict['score'])