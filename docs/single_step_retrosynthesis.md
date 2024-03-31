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
Aftet applying templates, this expander will search the original dataset by molecular similarity, providing reaction precedents.  
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
Expanders follow a standard interface, and can be imported and initialised as follows:  

(Please note, initialising an expander for the first time will automatically downloaded additional required files)  
```python
from rbc2 import RetroBioCatExpander, EnzymeMapExpander, BKMSExpander, RetroRulesExpander, AizynthFinderExpander, RingBreakerExpander, AskCosExpander

retrobiocat = RetroBioCatExpander()
enzymemap = EnzymeMapExpander()
bkms = BKMSExpander()
retrorules = RetroRulesExpander()
aizynthfinder = AizynthFinderExpander()
ringbreaker = RingBreakerExpander()
askcos = AskCosExpander()
```

## RetroBioCat expander - example usage

For stand-alone usage, all expanders have a .get_reactions(smi) method which takes a target SMILES does a single 
retrosynthesis step, producing [Reactions](data_model.md#reaction).

```python
from rbc2 import RetroBioCatExpander

expander = RetroBioCatExpander()
reactions = expander.get_reactions('CCCC=O')

# print the reaction smiles for the proposed reactions
for rxn in reactions:
    print(rxn.reaction_smiles())
    # CCCCO>>CCCC=O
    # CCCC(=O)O>>CCCC=O
    # CC=CC=O>>CCCC=O 
    # ect..

# Print all the precedents associated with the proposed reactions
for rxn in reactions:
    for precedent in rxn.precedents:
        print(precedent.data)  
```
See the [Data model](data_model.md) for further information on the attributes and methods which can be accessed for the proposed reactions.  

