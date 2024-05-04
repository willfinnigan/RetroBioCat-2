Using expanders as standalone tools
===================================


RetroBioCat expander
--------------------

::

    from rbc2 import RetroBioCatExpander

    expander = RetroBioCatExpander()
    reactions = expander.get_reactions('CCCC=O')

    for rxn in reactions:
        print(rxn.reaction_smiles())

outputs..
::

    CCCCO>>CCCC=O
    CCCC(=O)O>>CCCC=O
    CC=CC=O>>CCCC=O
    CCCCN>>CCCC=O
    CCCCN>>CCCC=O
    CCCC(O)C#N>>CCCC=O
    CCCC(=O)C(=O)O>>CCCC=O
    CCCC=N>>CCCC=O
