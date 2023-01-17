# Handle the workflow register
# This register contains processing details and warnings
# It may be useful to trace back some topology corrections
# WARNING: It is also used by some internal processes

import os
import json
from datetime import datetime

# Start the register dictionary with a few metadata
def start_register (inputs : dict) -> dict:
    return {
        'date': datetime.today().strftime('%d-%m-%Y %H:%M:%S'),
        'inputs': inputs,
        'warnings': [],
    }

# Set the register data filename
register_filename = 'register.json'

# Save the current register to a file
# In case there is a previous register, read it and append data to it
def save_register (register : 'Dependency'):
    register_data = []
    # Check if there is previous register data
    if os.path.exists(register_filename):
        with open(register_filename, 'r') as file:
            register_data = json.load(file)
    # Add current register data
    register_data.append(register)
    # Write it to a json file
    with open(register_filename, 'w') as file:
        json.dump(register_data, file, indent=4)