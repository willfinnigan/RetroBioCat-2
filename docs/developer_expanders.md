# Expanders

## Expander interface
All expanders in RetroBioCat 2.0 follow a standard interface, with the following methods:

```python
class Expander(ABC):
    """The expander interface, which defines the methods that an expander must implement"""

    @abstractmethod
    def __init__(self,
                 network: Optional[Network] = None,
                 config: Optional[Expansion_Config] = None):
        self.network = network
        self.rxn_type = ''  # eg retrobiocat, enzymemap, aizynthfinder ect..
        self.rxn_domain = ''  # biocatalysis, biosynthesis, or chemistry
        self.config = config

    @abstractmethod
    def get_reactions(self, smi: str) -> List[Reaction]:
        pass

    @abstractmethod
    def get_options(self, smi: str) -> List[ReactionOption]:
        pass

    @abstractmethod
    def create_option(self, smi: str, name: str, smarts: List[str],
                      template_metadata: dict, score: float) -> ReactionOption:
        pass

    @abstractmethod
    def number_of_rule_applications(self) -> int:
        pass

    @abstractmethod
    def number_of_calls(self) -> int:
        pass

```