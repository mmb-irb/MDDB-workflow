import json

from model_workflow.utils.file import File
from model_workflow.utils.structures import Structure
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
    if value_type == int or value_type == float: return value
    # For strings, sum the ordinal numbers of every letter
    if value_type == str: return sum(map(ord, value))
    # For objects, stringify them and then do the same that with strings
    if value_type in { list, dict }:
        stringifyed = json.dumps(value)
        return get_cksum_id(stringifyed)
    # For files use file last modification time and size
    if isinstance(value, File): return f'{value.mtime}-{value.size}'
    # For the parsed structure
    if isinstance(value, Structure):
        pdb_content = value.generate_pdb()
        return get_cksum_id(pdb_content)
    # If the value has non of previous types then we complain
    raise TypeError(f'Non supported type "{value_type}" for cksum id: {value}')