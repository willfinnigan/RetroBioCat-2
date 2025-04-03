import json
from pathlib import Path

from rbc2.data_model.reaction import reactions_to_dicts
from rbc2.pathway_tools.pa_route_conversion import reactions_from_pa_route

test_data = Path(__file__).parents[1] / 'test_data'

def test_load_pathway_paroutes():
    filepath = f"{test_data}/pa_route.json"
    with open(filepath, 'r') as f:
        route_data = json.load(f)

    print(type(route_data))

    print(route_data)

    reactions = reactions_from_pa_route(route_data)

    print(reactions)

    for rxn in reactions:
        print(rxn)

    reaction_dicts = reactions_to_dicts(reactions)