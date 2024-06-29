# The data model

RetroBioCat 2.0 uses the following dataclasses to hold the results of automated retrosynthesis:

## Reaction

The `Reaction` class is used to represent a single reaction.  

#### Attributes
`product: str` - The SMILES for the product of the reaction  
`substrates: List[str]` - A list of SMILES for the substrates of the reaction  
`score: float` - A score given to the reaction by the expander which created it.  Expanders have differnt scoring methods.  
`name: str` - A name given to the reaction, usually by the expansion method which created it  
`rxn_type: str` - The type of expander which created this reaction (eg, 'retrobiocat', 'aizynthfinder' ect)  
`rxn_domain: str` - The domain of the reaction - currently either 'biocatalysis', 'chemistry' or 'biosynthesis'  
`unique_id: str` - A unique identifier for the reaction, normally uuid4  
`template_metadata: dict[str: dict]` - A dictionary of metadata about the template used to create the reaction  
`precedents: List[Precedent]` - A list of Precedent objects, which are similar reactions in the literature  
`feasability_filter_scores: dict[str: float]` - A dictionary of scores from different feasability filters  

#### Methods
`get_complexity_change() -> float` - Returns the change in complexity of the reaction  
`reaction_smiles() -> str` - Returns a SMILES string representation of the reaction.
`get_similarity_score() -> float` - Returns a similarity score for the reaction, from the closest precedent  
`to_dict() -> dict` - Returns a dictionary representation of the reaction  

#### Other standalone functions


---

## Precedent
Reactions can have precedents, which are typically similar reactions that have been seen before.  
Each reaction has a `precedents` attribute which is a list of `Precedent` objects.

#### Attributes
`name: str` - The name of the precedent  
`rxn_smi: str` - The SMILES string of the precedent reaction  
`precedent_id: str` - A unique identifier for the precedent  
`data: dict` - A dictionary of metadata about the precedent - depends on the source of the precedent  
`similarity: float` - A similarity score between the precedent and the reaction - typically on the propducts    

---

## Pathway
A `Pathway` is a list of `Reaction` objects, representing a series of reactions to get from a starting material to a target product.

```python
from rbc2 import Pathway

pathway = Pathway([reaction1, 
                   reaction2, 
                   reaction3])
```

#### Attributes
`reactions: List[Reaction]` - A list of `Reaction` objects  
`smi_produced_by: dict[str: Set[Reaction]]` - For every SMILES in the pathway, access the reaction(s) which produce it  
`smi_substrate_of: dict[str: Set[Reaction]]` - For every SMILES in the pathway, access the reaction(s) which use it as a substrate  
`product_smis: Set[str]` - A set of all the SMILES which are products of any reaction in the pathway  
`substrate_smis: Set[str]` - A set of all the SMILES which are substrates of any reaction in the pathway  
`all_smis: Set[str]` - A set of all the SMILES in the pathway  
`target_smi: str` - The target SMILES of the pathway - this is the product of the final reaction, can be passed in as a parameter on init.  
`pathway_length: int` - How many steps deep does the pathway go.  != len(reactions) if there are branches.  
`end_smi_depths: dict[str: int]` - A dictionary of the depth of each SMILES in the pathway  
`tree: dict` - A tree representation of the pathway, starting from target_smi in the format {'smiles': smi, 'depth': depth, 'children': []}  

#### Methods
`get_pa_route(starting_material_evaluator) -> dict` - Output the pathway in the format used by PARoutes. 
Needs a [starting_material_evaluator](starting_materials.md) to be passed in to evaluate what molecules are available.  
`end_smis() -> List[str]` - Return a list of the SMILES which are substrates of the initial reaction(s) in the pathway  
`save() -> List[dict]` - Save the pathway to a json file.  Returns the reactions in the pathway as a list of dictionaries.  
`get_reaction_with_product(smi: str) -> Reaction:` - Return the first reaction which produces the SMILES.  
`get_reaction_with_substrate(smi: str) -> Reaction:` - Return the first reaction which uses the SMILES as a substrate.  

#### Other standalone functions
`load_pathway(reaction_dict_list: List[dict]) -> Pathway` - Load a pathway from a list of dictionaries.  Loads the format which is saved by pathway.save()    




