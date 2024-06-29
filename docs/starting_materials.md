# Starting Material Evaluator


## Commercial Starting Material Evaluator (default)
RetroBioCat 2 has a default starting material evaluator, which utilises a database of commercially available molecules.

Initialising this evaluator for the first time will download the necessary sqlite database file.

```python
from rbc2 import CommercialSME

sme = CommercialSME()
available, info = sme.eval('CC(=O)Oc1ccccc1C(=O)O')
print(available, info)
```

## E coli metabolism Starting Material Evaluator
RetroBioCat 2 also has a starting material evaluator for E coli metabolism.

```python
from rbc2 import EcoliSME

sme = EcoliSME()
available, info = sme.eval('CC(=O)Oc1ccccc1C(=O)O')
print(available, info)
```



## Starting Material Evaluator Interface

An abstract base class exists for starting material evaluators, which defines the expected methods.

A custom evaluator can therefore be created by subclassing this interface, and utilised in the same way as the default evaluator.

```
class StartingMaterialEvaluatorInterface(ABC):

    def __init__(self):
        self.target_always_not_buyable = True

    @abstractmethod
    def eval(self, smi: str) -> Tuple[bool, dict]:
        """Evaulate the smiles for whether its available as a starting material.
        Returns a tuple of (bool, dict) where the bool is True if the smiles is available and the dict contains
        additional information about the availability.
        """
        pass

    @abstractmethod
    def is_mol_chiral(self, smi: str) -> bool:
        """Return True if the molecule is chiral"""
        pass
```
