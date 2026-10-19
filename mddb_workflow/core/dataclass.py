from dataclasses import dataclass, field
from mddb_workflow.utils.type_hints import *


@dataclass
class InChIKeyData:
    """Data structure for InChI key information.

    Attributes:
        inchi (str): The InChI string for the residue.
        resindices (list[int]): List of all residue indices with this InChI key.
        fragments (list[list[int]]): Lists of residue indices that are connected as a single group.
        resnames (set[str]): Set of residue names associated with this InChI key.
        molname (str): Representative molecule name for this InChI key.
        moltype (Literal['residue', 'fragment']): Type of the molecule.
        classification (set): Set of residue classifications for this InChI key.
        frag_len (int): Length on number of residues of the fragments. 1 if no fragments are present.
        references (dict): Additional database references related to this InChI key.

    """
    inchi: str
    resindices: list[int] = field(default_factory=list)
    fragments: list[list[int]] = field(default_factory=list)
    resnames: set[str] = field(default_factory=set)
    molname: str = ''
    moltype: Literal['residue', 'fragment'] = 'residue'
    classification: set = field(default_factory=set)
    frag_len: int = 1
    references: dict = field(default_factory=dict)

    @classmethod
    def load_cache(cls, cache_dict: dict[dict | None]) -> dict[str, 'InChIKeyData']:
        """Load data from cache by converting dictionaries to InChIKeyData objects."""
        if len(cache_dict) > 0 and isinstance(next(iter(cache_dict.values())), dict):
            inchikeys = {}
            for inchikey, data in cache_dict.items():
                inchidata = cls(inchi=data['inchi'])
                for key, value in data.items():
                    setattr(inchidata, key, value)
                inchikeys[inchikey] = inchidata
            return inchikeys
        else:
            return cache_dict
