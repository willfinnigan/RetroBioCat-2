from rbc2 import RetroBioCatExpander

if __name__ == '__main__':
    expander = RetroBioCatExpander(include_experimental=True,
                                   include_two_step=True,
                                   include_requires_absence_of_water=True)
    expander.rxn_class.get_rxns()
    expander.rxn_class.get_multistep_rxns()
    rxn_enzyme_map = expander.rxn_class.reaction_enzyme_map

    print(rxn_enzyme_map)

    # invert the rxn_enzyme_map so its is enzyme -> rxns
    enzyme_rxn_map = {}
    for rxn, enzymes in rxn_enzyme_map.items():
        for enzyme in enzymes:
            if enzyme not in enzyme_rxn_map:
                enzyme_rxn_map[enzyme] = []
            enzyme_rxn_map[enzyme].append(rxn)

    print(enzyme_rxn_map)

    # save as a json
    import json
    with open('enzyme_rxn_map.json', 'w') as f:
        json.dump(enzyme_rxn_map, f, indent=4)




