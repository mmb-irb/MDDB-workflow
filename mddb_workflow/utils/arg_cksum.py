import json

from mddb_workflow.utils.auxiliar import safe_hasattr
from mddb_workflow.utils.file import File
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.register import Register
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.mda_spells import get_mda_universe_cksum
from mddb_workflow.utils.type_hints import *

from MDAnalysis.core.universe import Universe

# Note that it is more convinient not to have this function among auxiliar functions
# We depend on class types which depend on the auxiliar module

# Given a value of any type, generate a 'cksum' like code which is reproducible
# This is used to compare values between different runs without having to store the whole value
def get_cksum_id (value) -> Optional[int | float | str]:
    # Nones remain as they are
    if value == None: return None
    # Ge the value type
    value_type = type(value)
    # For numbers simply use the number itself
    if value_type == int or value_type == float or value_type == bool: return value
    # For strings, sum the ordinal numbers of every letter
    if value_type == str: return f'{sum(map(ord, value))} -> {len(value)}'
    # For exceptions simply keep the excepction message
    if value_type == Exception: return str(value)
    # For objects, stringify them and then do the same that with strings
    if value_type in { list, dict, tuple }:
        stringifyed = json.dumps(value, default = lambda o: o.__repr__())
        return get_cksum_id(stringifyed)
    # If we have a set then make a list and sort it alphabetically
    # Thus we make sure the order is coherent between different
    if value_type == set:
        standard = list(value).sort()
        return get_cksum_id(standard)
    # For functions simply store the name
    if callable(value): return value.__name__
    # For files use file last modification time and size
    if isinstance(value, File): return value.get_cksum()
    # For the parsed structure
    if isinstance(value, Structure):
        pdb_content = value.generate_pdb()
        return get_cksum_id(pdb_content)
    # For a MDAnalysis universe
    if isinstance(value, Universe): return get_mda_universe_cksum(value)
    # For the parsed structure
    if isinstance(value, Selection): return f'{len(value.atom_indices)}-{sum(value.atom_indices)}'
    # For handler class instances it makes not sense making a comparision
    # WARNING: If they are used as input then make sure it is only being "used" and not "read"
    # Otherwise, the content will be not compared between runs (e.g. warnings in the register)
    if isinstance(value, Register): return True
    if isinstance(value, Cache): return True
    if safe_hasattr(value, '__class__'):
        if value.__class__.__name__ == 'Project': return True
        if value.__class__.__name__ == 'MD': return True
    # If the value has non of previous types then we complain
    raise TypeError(f'Non supported type "{value_type}" for cksum id: {value}')