from sys import argv

from mddb_wf.utils.auxiliar import load_yaml, save_yaml

# Read a reference inputs file
reference_inputs_path = 'inputs.yaml'
reference_inputs = load_yaml(reference_inputs_path)

# Get the directories where this logic is to be run
if len(argv) < 2:
    raise SystemExit('Missing target directories')

# Iterate target directories
for directory in argv[1:]:
    # Set custom input values for this directory
    # Copy the reference inputs
    new_inputs = { k: v for k, v in reference_inputs.items() }
    # Update the pdb ids and the title
    pdb_id = directory.replace('/','').replace('\n','')
    new_inputs['pdb_ids'] = [pdb_id]
    new_inputs['name'] += ', ' + pdb_id
    new_inputs_path = directory + 'inputs.yaml'
    save_yaml(new_inputs, new_inputs_path)