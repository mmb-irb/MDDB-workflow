import json

from model_workflow.utils.file import File
from model_workflow.utils.selections import Selection
from model_workflow.utils.structures import Structure
from model_workflow.utils.register import Register
from model_workflow.utils.cache import Cache
from typing import Optional, Union

# Note that it is more convinient not to have this function among auxiliar functions
# We depend on class types which depend on the auxiliar module

# Given a value of any type, generate a 'cksum' like code which is reproducible
# This is used to compare values between different runs without having to store the whole value
def get_cksum_id (value) -> Optional[Union[int, float, str]]:
    # Nones remain as they are
    if value == None: return None
    # Ge the value type
    value_type = type(value)
    # For numbers simply use the number itself
    if value_type == int or value_type == float or value_type == bool: return value
    # For strings, sum the ordinal numbers of every letter
    if value_type == str: return f'{sum(map(ord, value))} -> {len(value)}'
    # For objects, stringify them and then do the same that with strings
    if value_type in { list, dict }:
        stringifyed = json.dumps(value, default = lambda o: '<not serializable>')
        return get_cksum_id(stringifyed)
    # For functions simply store the name
    if callable(value): return value.__name__
    # For files use file last modification time and size
    if isinstance(value, File):
        if not value.exists: return f'missing {value.path}'
        return f'{value.path} -> {value.mtime} {value.size}'
    # For the parsed structure
    if isinstance(value, Structure):
        pdb_content = value.generate_pdb()
        return get_cksum_id(pdb_content)
    # For the parsed structure
    if isinstance(value, Selection): return f'{len(value.atom_indices)}-{sum(value.atom_indices)}'
    # For the register or the cache it makes not sense making a comparision
    # WARNING: If the register/cache is an input then make sure it is only being "used" and not "read"
    # Otherwise, the register/cache content (e.g. warnings) will be not compared between runs
    # If you want this to be compared then pass the register/cache sub-value as an input
    if isinstance(value, Register): return True
    if isinstance(value, Cache): return True
    # If the value has non of previous types then we complain
    raise TypeError(f'Non supported type "{value_type}" for cksum id: {value}')