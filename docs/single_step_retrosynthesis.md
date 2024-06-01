# Expanders for single step retrosynthesis

Expanders are objects which perform single step retrosynthesis - 
that is, given a target molecule they will return a list of reactions which can be used to make that molecule. 
The reactions are represented as [Reaction](data_model.md#reaction) objects.

RetroBioCat 2.0 has a range of expanders available, including of course the RetroBioCat expander, 
which applies manually curated reaction rules for the biocatalysis toolbox.

## Enzyme expanders
**RetroBioCat**  
Manually curated reaction rules which capture the enzyme toolbox for biocatlysis.  Searches the retrobiocat database (not included by default) for reaction precedents.  
`Finnigan, W., Hepworth, L. J., Flitsch, S. L. & Turner, N. J. RetroBioCat as a computer-aided synthesis planning tool for biocatalytic reactions and cascades. Nature Catalysis 4, 98–104 (2021).`  

**EnzymeMap**  
EnzymeMap is a set of scripts to process and clean BRENDA, followed by automated template extraction and the creation of a template relevance model. 
After applying templates, this expander will search the original dataset by molecular similarity, providing reaction precedents.   
`Heid, E., Probst, D., Green, W. H. & Madsen, G. K. H. EnzymeMap: curation, validation and data-driven prediction of enzymatic reactions. Chem. Sci. 14, 14229–14242 (2023).`  

**BKMS**  
Data from the BKMS database processed into SMILES format, followed by automated template extraction and the creation of a template relevance model. 
After applying templates, this expander will search the original dataset by molecular similarity, providing reaction precedents.  
`Levin, I., Liu, M., Voigt, C. A. & Coley, C. W. Merging enzymatic and synthetic chemistry with computational synthesis planning. Nat Commun 13, 7747 (2022).`  

**RetroRules**  
RetroRules templates used by RetroPath.  These templates have been automatically extracted at a range of diameters. 
This expander prioritises templates by scoring them using a combination of product similarity and a predetermined biological score.  
  `Koch, M., Duigou, T. & Faulon, J.-L. Reinforcement Learning for Bioretrosynthesis. ACS Synth. Biol. 9, 157–168 (2020).`  

## Chemistry expanders
**AIZynthfinder**  
`Genheden, S. et al. AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. J Cheminform 12, 70 (2020).`  

**RingBreaker**  
`Thakkar, A., Selmi, N., Reymond, J.-L., Engkvist, O. & Bjerrum, E. J. “Ring Breaker”: Neural Network Driven Synthesis Prediction of the Ring System Chemical Space. J. Med. Chem. 63, 8791–8808 (2020).`  

**AskCos**  
`Coley, C. W. et al. A robotic platform for flow synthesis of organic compounds informed by AI planning. Science 365, eaax1566 (2019).`

## Usage
Expanders follow a standard interface, and can be imported and initialised as follows. 

Default keyword arguments are shown. These can be excluded to use the default values, which is often the best option.  

(Please note, initialising an expander for the first time will automatically downloaded additional required files)  
```python
from rbc2 import RetroBioCatExpander, EnzymeMapExpander, BKMSExpander, RetroRulesExpander, AIZynthfinderExpander, RingBreakerPolicyExpander, AskcosPolicyExpander

retrobiocat = RetroBioCatExpander(include_experimental = False,
                                  include_two_step = True,
                                  include_requires_absence_of_water = False,
                                  score_similarity_before_option_creation  = True,
                                  search_literature_precedent  = True,
                                  only_active_literature_precedent  = True,
                                  similarity_cutoff = 0.55)

enzymemap = EnzymeMapExpander(cutoff_cumulative=0.995,
                              cutoff_number=50,
                              enable_precedent_search=True,
                              similarity_cutoff=0.1)


bkms = BKMSExpander(cutoff_cumulative=0.995,
                    cutoff_number=50,
                    allow_multi_product_templates=False,
                    enable_precedent_search=True,
                    similarity_cutoff=0.1)

retrorules = RetroRulesExpander(rank_by='combined_score',  # options are: similarity, score, combined_score
                                diameter=6,
                                similarity_threshold=0.2,
                                score_threshold=0.2,
                                combined_score_threshold=0.2,
                                max_reactions=100) # max number of similar substrates to consider reactions for)


aizynthfinder = AIZynthfinderExpander(cutoff_cumulative=0.995,
                                      cutoff_number=50)

ringbreaker = RingBreakerPolicyExpander(cutoff_cumulative=0.995,
                                  cutoff_number=10)

askcos = AskcosPolicyExpander(cutoff_cumulative=0.995,
                        cutoff_number=50)
```

## RetroBioCat expander - example usage

For stand-alone usage, all expanders have a .get_reactions(smi) method which takes a target SMILES and does a single 
retrosynthesis step, producing [Reactions](data_model.md#reaction).

```python
from rbc2 import RetroBioCatExpander

expander = RetroBioCatExpander()
reactions = expander.get_reactions('CCCC=O')

# print the reaction smiles for the proposed reactions
for rxn in reactions:
    print(rxn.reaction_smiles())  # eg CCCCO>>CCCC=O

# Print all the precedents associated with the proposed reactions
for rxn in reactions:
    for precedent in rxn.precedents:
        print(precedent.data)  
```
See the [Data model](data_model.md) for further information on the attributes and methods which can be accessed for the proposed reactions.  

## Expander Configuration
As well as the specific keyword arguments available for each expander, 
all expanders have a `config` argument which can be used to pass in a configuration object for expanders in general.  

```python
from rbc2.configs import Expansion_Config
from rbc2 import RetroBioCatExpander

config = Expansion_Config()
config.max_reactions = 10   # for example, change the maximum number of reactions returned by any expander

expander = RetroBioCatExpander(config=config)
```

See all the options and their defaults at [source](https://github.com/willfinnigan/RetroBioCat_2/blob/main/rbc2/configs/expansion_config.py)