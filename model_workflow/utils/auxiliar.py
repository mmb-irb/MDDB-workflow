# Auxiliar generic functions and classes used along the workflow

from model_workflow.utils.constants import RESIDUE_NAME_LETTERS

import sys

# Check if a module has been imported
def is_imported (module_name : str) -> bool:
    return module_name in sys.modules

# Set a custom exception for user input errors to avoid showing traceback in the terminal
class QuietException (Exception):
    pass

# Set custom exception which are to be quiet because they are not caused by an error in the code
class InputError (QuietException):
    pass

class TestFailure (QuietException):
    pass

class MissingDependency (QuietException):
    pass

# Set a custom exception handler where our input error exception has a quiet behaviour
def custom_excepthook (exception_class, message, traceback):
    # Quite behaviour if it is our input error exception
    if QuietException in exception_class.__bases__:
        print('{0}: {1}'.format(exception_class.__name__, message))  # Only print Error Type and Message
        return
    # Default behaviour otherwise
    sys.__excepthook__(exception_class, message, traceback)
sys.excepthook = custom_excepthook

# Set a function to get the next letter from an input letter in alphabetic order
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
def get_new_letter(current_letters : list) -> str:
    new_letter = next((letter for letter in letters if letter not in current_letters), None)
    if not new_letter:
        raise Exception("There are no more letters in the alphabet")
    return new_letter

# Given a residue name, return its single letter
def residue_name_to_letter (residue_name : str) -> str:
    return RESIDUE_NAME_LETTERS.get(residue_name, 'X')