# Multi-step retrosynthesis using MCTS

For planning multi-step routes, RetroBioCat 2 uses a Monte Carlo Tree Search (MCTS) algorithm.

Minimal quick start example
```python
from rbc2 import MCTS
from rbc2 import get_expanders

 # Target SMILES - this is what we want to make
target_smi = 'CCCC=O' 

# Get the expanders with a convenience function
expanders = get_expanders(['retrobiocat', 
                           'aizynthfinder'])  

# Initialise and run the MCTS
mcts = MCTS(target_smi, expanders)
mcts.config.max_search_time = 15  # only search for 15 seconds
mcts.run()  

# Get the output. This will be a list of Pathways
all_pathways = mcts.get_all_pathways()
solved_pathways = mcts.get_solved_pathways()
```





