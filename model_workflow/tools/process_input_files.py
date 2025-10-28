from model_workflow.utils.auxiliar import InputError, MISSING_TOPOLOGY, warn
from model_workflow.utils.constants import STRUCTURE_FILENAME, TRAJECTORY_FILENAME
from model_workflow.utils.constants import CONVERTED_STRUCTURE, CONVERTED_TRAJECTORY
from model_workflow.utils.constants import FILTERED, FILTERED_STRUCTURE, FILTERED_TRAJECTORY
from model_workflow.utils.constants import IMAGED, IMAGED_STRUCTURE, IMAGED_TRAJECTORY
from model_workflow.utils.constants import CORRECTED_STRUCTURE, CORRECTED_TRAJECTORY
from model_workflow.utils.constants import INCOMPLETE_PREFIX, CG_ATOM_ELEMENT, SNAPSHOTS_FLAG
from model_workflow.utils.file import File
from model_workflow.utils.structures import Structure
from model_workflow.utils.pyt_spells import get_frames_count
from model_workflow.utils.arg_cksum import get_cksum_id
from model_workflow.utils.type_hints import *

from model_workflow.tools.check_inputs import check_inputs, PREFILTERED_TOPOLOGY_EXCEPTION
from model_workflow.tools.conversions import convert
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.structure_corrector import structure_corrector


def is_amber_top (input_topology_file : 'File') -> bool:
    """ Check if a .top file is from Amber. 
    Returns True if it is Amber, False if it is Gromacs. """
    if input_topology_file != MISSING_TOPOLOGY and \
        input_topology_file.extension == 'top':
        with open(input_topology_file.path, 'r') as f:
            lines = f.readlines(5)

            # If all non-empty first words are '%' assume Amber (.prmtop)
            first_words = {line.split()[0] for line in lines if line.strip()}
            if '%VERSION' in first_words:
                return True

            # If any line starts with ';' or '[' assume Gromacs (.top)
            first_chars = {word[0] for word in first_words}
            if any(c in (';', '[') for c in first_chars):
                return False

            # Otherwise we cannot decide
            raise InputError('Unable to infer topology format from first five lines')
            
    return False


def process_input_files (
    input_structure_file : 'File',
    input_trajectory_files : 'File',
    input_topology_file : 'File',
    # Output
    output_directory : str,
    topology_filepath : str,
    # Processing parameters
    filter_selection : str,
    image : bool,
    fit : bool,
    translation : tuple,
    # Make sure the MD is used only to set values or use its functions, but not to get values
    # Values msut be passed separatedly as inputs so the taks can identify when inputs change
    self : 'MD',
    # Get the task which is calling this function
    # Thus we may knwo which inputs have changed compared to previous runs
    task : 'Task',
):
    """ Process input files to generate the processed files.
    This process corrects and standarizes the topology, the trajectory and the structure. """

    # Make sure we do not enter in a loop
    # This may happen when we read/call an output value/file by mistake
    if hasattr(self, '_processed'): raise RuntimeError('Looped processing')
    self._processed = True

    # Input trajectories should have all the same format
    input_trajectory_formats = set([ trajectory_file.format for trajectory_file in input_trajectory_files ])
    if len(input_trajectory_formats) > 1:
        raise InputError('All input trajectory files must have the same format')

    # Set the output filepaths
    # Note that these outputs aim for the MD directory, not for the task directory
    output_structure_filepath = self.pathify(STRUCTURE_FILENAME)
    output_structure_file = File(output_structure_filepath)
    output_trajectory_filepath = self.pathify(TRAJECTORY_FILENAME)
    output_trajectory_file = File(output_trajectory_filepath)
    output_topology_filepath = topology_filepath
    output_topology_file = File(output_topology_filepath) if output_topology_filepath else MISSING_TOPOLOGY

    # --- TOPOLOGY FORMAT ASSUMTION ---------------------------------------------------------

    # Make a super fast check and an assumption
    # Topologies with the .top extension for us are considered Gromacs topology format
    # However it is usual than Amber topologies (ideally .prmtop) are also '.top'
    # So if the trajectory is Amber and the topology is .top then assume it is Amber
    input_trajectories_format = list(input_trajectory_formats)[0]
    if is_amber_top(input_topology_file):
        # Creating a topology symlink/copy with the correct extension
        warn(f'Topology is .top but the trajectory is from Amber. It is assumed the topology is .prmtop')
        reformatted_topology_file = input_topology_file.reformat('prmtop')
        output_topology_file.path = output_topology_file.extensionless_filepath+'.prmtop'
        if input_structure_file == input_topology_file:
            input_structure_file = reformatted_topology_file
        input_topology_file = reformatted_topology_file
        
    # --- FIRST CHECK -----------------------------------------------------------------------

    # Check input files match in number of atoms
    # Here we have not standarized the format so we must check differently with every format
    exceptions = check_inputs(input_structure_file, input_trajectory_files, input_topology_file)

    # There is a chance that the inputs checker has prefiltered the topology to match trajectory
    # If this is the case then use the prefiltered topology from now on
    prefiltered_topology = exceptions.get(PREFILTERED_TOPOLOGY_EXCEPTION, None)
    if prefiltered_topology:
        if input_structure_file == input_topology_file:
            input_structure_file = prefiltered_topology
        input_topology_file = prefiltered_topology

    # --- CONVERTING AND MERGING ------------------------------------------------------------

    # Set the output format for the already converted structure
    input_structure_format = input_structure_file.format
    output_structure_format = output_structure_file.format
    converted_structure_filepath = f'{output_directory}/{CONVERTED_STRUCTURE}'
    # If input structure already matches the output format then avoid the renaming
    if input_structure_format == output_structure_format:
        converted_structure_filepath = input_structure_file.path
    # Set the output file for the already converted structure
    converted_structure_file = File(converted_structure_filepath)
    # Set the output format for the already converted trajectory
    input_trajectories_format = list(input_trajectory_formats)[0]
    output_trajectory_format = output_trajectory_file.format
    # Set the output file for the already converted trajectory
    converted_trajectory_filepath = f'{output_directory}/{CONVERTED_TRAJECTORY}'
    # If input trajectory already matches the output format and is unique then avoid the renaming
    if input_trajectories_format == output_trajectory_format and len(input_trajectory_files) == 1:
        converted_trajectory_filepath = input_trajectory_files[0].path
    converted_trajectory_file = File(converted_trajectory_filepath)
    # Join all input trajectory paths
    input_trajectory_paths = [ trajectory_file.path for trajectory_file in input_trajectory_files ]

    # Set an intermeidate file for the trajectory while it is being converted
    # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while converting
    incompleted_converted_trajectory_filepath = f'{output_directory}/{INCOMPLETE_PREFIX + CONVERTED_TRAJECTORY}'
    incompleted_converted_trajectory_file = File(incompleted_converted_trajectory_filepath)
    # If there is an incomplete trajectory then remove it
    if incompleted_converted_trajectory_file.exists:
        incompleted_converted_trajectory_file.remove()

    # Convert input structure and trajectories to output structure and trajectory
    if not converted_structure_file.exists or not converted_trajectory_file.exists:
        print(' * Converting and merging')
        convert(
            input_structure_filepath = input_structure_file.path,
            output_structure_filepath = converted_structure_file.path,
            input_trajectory_filepaths = input_trajectory_paths,
            output_trajectory_filepath = incompleted_converted_trajectory_file.path,
        )
        # Once converted, rename the trajectory file as completed
        # If the converted trajectory already exists then it means it is the input trajectory
        if converted_trajectory_file.exists:
            incompleted_converted_trajectory_file.remove()
        else:
            incompleted_converted_trajectory_file.rename_to(converted_trajectory_file)

    # Topologies are never converted, but they are kept in their original format

    # --- provisional reference structure ---

    # Now that we MUST have a PDB file we can set a provisional structure instance
    # Note that this structure is not yet corrected so it must be used with care
    # Otherwise we could have silent errors
    provisional_structure = Structure.from_pdb_file(converted_structure_file.path)
    # Now we can set a provisional coarse grain selection
    # This selection is useful to avoid problems with CG atom elements
    # Since this is proviosonal we will make it silent
    provisional_cg_selection = self._set_cg_selection(provisional_structure, verbose=False)
    for atom_index in provisional_cg_selection.atom_indices:
        provisional_structure.atoms[atom_index].element = CG_ATOM_ELEMENT

    # --- FILTERING ATOMS ------------------------------------------------------------

    # Find out if we need to filter
    # i.e. check if there is a selection filter and it matches some atoms
    must_filter = bool(filter_selection)

    # Set output filenames for the already filtered structure and trajectory
    # Note that this is the only step affecting topology and thus here we output the definitive topology
    filtered_structure_filepath = f'{output_directory}/{FILTERED_STRUCTURE}'
    if not must_filter:
        filtered_structure_filepath = converted_structure_filepath
    filtered_structure_file = File(filtered_structure_filepath)
    filtered_trajectory_filepath = f'{output_directory}/{FILTERED_TRAJECTORY}'
    if not must_filter:
        filtered_trajectory_filepath = converted_trajectory_filepath
    filtered_trajectory_file = File(filtered_trajectory_filepath)
    filtered_topology_file = output_topology_file if must_filter else input_topology_file

    # Set an intermeidate file for the trajectory while it is being filtered
    # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while filtering
    incompleted_filtered_trajectory_filepath = f'{output_directory}/{INCOMPLETE_PREFIX + FILTERED_TRAJECTORY}'
    incompleted_filtered_trajectory_file = File(incompleted_filtered_trajectory_filepath)
    # If there is an incomplete trajectory then remove it
    if incompleted_filtered_trajectory_file.exists:
        incompleted_filtered_trajectory_file.remove()

    # Check if any output file is missing
    missing_filter_output = not filtered_structure_file.exists \
        or not filtered_trajectory_file.exists \
        or (filtered_topology_file != MISSING_TOPOLOGY and not filtered_topology_file.exists)

    # Check if parameters have changed
    # Note that for this specific step only filtering is important
    previous_filtered_parameters = self.cache.retrieve(FILTERED)
    same_filtered_parameters = previous_filtered_parameters == filter_selection
    
    # Filter atoms in structure, trajectory and topology if required and not done yet
    if must_filter and (missing_filter_output or not same_filtered_parameters):
        print(' * Filtering atoms')
        filter_atoms(
            input_structure_file = converted_structure_file,
            input_trajectory_file = converted_trajectory_file,
            input_topology_file = input_topology_file, # We use input topology
            output_structure_file = filtered_structure_file,
            output_trajectory_file = incompleted_filtered_trajectory_file,
            output_topology_file = filtered_topology_file, # We genereate the definitive topology
            reference_structure = provisional_structure,
            filter_selection = filter_selection,
        )
        # Once filetered, rename the trajectory file as completed
        # If the filtered trajectory already exists then it means it is the input trajectory
        if filtered_trajectory_file.exists:
            incompleted_filtered_trajectory_file.remove()
        else:
            incompleted_filtered_trajectory_file.rename_to(filtered_trajectory_file)
        # Update the cache
        self.cache.update(FILTERED, filter_selection)

    # --- provisional reference structure ---

    # Now that we have a filtered PDB file we have to update provisional structure instance
    # Note that this structure is not yet corrected so it must be used with care
    # Otherwise we could have silent errors
    provisional_structure = Structure.from_pdb_file(filtered_structure_file.path)
    # Again, set the coarse grain atoms
    # Since elements may be needed to guess PBC selection we must solve them right before
    # Since this is proviosonal we will make it silent
    provisional_cg_selection = self._set_cg_selection(provisional_structure, verbose=False)
    for atom_index in provisional_cg_selection.atom_indices:
        provisional_structure.atoms[atom_index].element = CG_ATOM_ELEMENT
    # Also we can set a provisional PBC selection
    # This selection is useful both for imaging/fitting and for the correction
    # We will make sure that the provisonal and the final PBC selections match
    # Since this is proviosonal we will make it silent
    provisional_pbc_selection = self._set_pbc_selection(provisional_structure, verbose=False)

    # --- IMAGING AND FITTING ------------------------------------------------------------

    # There is no logical way to know if the trajectory is already imaged or it must be imaged
    # We rely exclusively in input flags
    must_image = image or fit

    # Set output filenames for the already filtered structure and trajectory
    imaged_structure_filepath = f'{output_directory}/{IMAGED_STRUCTURE}'
    if not must_image:
        imaged_structure_filepath = filtered_structure_filepath
    imaged_structure_file = File(imaged_structure_filepath)
    imaged_trajectory_filepath = f'{output_directory}/{IMAGED_TRAJECTORY}'
    if not must_image:
        imaged_trajectory_filepath = filtered_trajectory_filepath
    imaged_trajectory_file = File(imaged_trajectory_filepath)

    # Set an intermeidate file for the trajectory while it is being imaged
    # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while imaging
    incompleted_imaged_trajectory_filepath = f'{output_directory}/{INCOMPLETE_PREFIX + IMAGED_TRAJECTORY}'
    incompleted_imaged_trajectory_file = File(incompleted_imaged_trajectory_filepath)
    # If there is an incomplete trajectory then remove it
    if incompleted_imaged_trajectory_file.exists:
        incompleted_imaged_trajectory_file.remove()

    # Check if any output file is missing
    missing_imaged_output = not imaged_structure_file.exists or not imaged_trajectory_file.exists

    # Check if parameters have changed
    # Note that for this step the filter parameters is also important
    previous_imaged_parameters = self.cache.retrieve(IMAGED)
    same_imaged_parameters = previous_imaged_parameters == (image, fit, *translation)
    print(f'IMAGE PARAMETERS: {previous_imaged_parameters} vs {(image, fit, *translation)} -> {same_imaged_parameters}')

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    if must_image and (missing_imaged_output or not same_imaged_parameters):
        print(' * Imaging and fitting')
        image_and_fit(
            input_structure_file = filtered_structure_file,
            input_trajectory_file = filtered_trajectory_file,
            input_topology_file = filtered_topology_file, # This is optional if there are no PBC residues
            output_structure_file = imaged_structure_file,
            output_trajectory_file = incompleted_imaged_trajectory_file,
            image = image,
            fit = fit,
            translation = translation,
            structure = provisional_structure,
            pbc_selection = provisional_pbc_selection
        )
        # Once imaged, rename the trajectory file as completed
        # If the imaged trajectory already exists then it means it is the input trajectory
        if imaged_trajectory_file.exists:
            incompleted_imaged_trajectory_file.remove()
        else:
            incompleted_imaged_trajectory_file.rename_to(imaged_trajectory_file)
        # Update the cache
        self.cache.update(IMAGED, (image, fit, *translation))
        # Update the provisional strucutre coordinates
        imaged_structure = Structure.from_pdb_file(imaged_structure_file.path)
        imaged_structure_coords = [ atom.coords for atom in imaged_structure.atoms ]
        provisional_structure.set_new_coordinates(imaged_structure_coords)

    # --- CORRECTING STRUCTURE ------------------------------------------------------------

    # Note that this step, although it is foucsed in the structure, requires also the trajectory
    # Also the trajectory may be altered in very rare cases where coordinates must be resorted

    # There is no possible reason to not correct the structure
    # This is the last step so the output files will be named as the output files of the whole processing

    # WARNING:
    # For the correcting function we need the number of snapshots and at this point it should not be defined
    # Snapshots are calculated by default from the already processed structure and trajectory
    # For this reason we can not rely on the public snapshots getter
    # We must calculate snapshots here using last step structure and trajectory
    snapshots = None
    # If we already have a value in the cache then use it, unless the input trajectory has changed
    same_trajectory = 'input_trajectory_files' not in task.changed_inputs
    if same_trajectory: snapshots = self.cache.retrieve(SNAPSHOTS_FLAG)
    # Calculate the new value
    if snapshots == None:
        snapshots = get_frames_count(imaged_structure_file, imaged_trajectory_file)
    # Update the MD snapshots value
    self.get_snapshots._set_parent_output(self, snapshots)
    # Save the snapshots value in the cache as well
    self.cache.update(SNAPSHOTS_FLAG, snapshots)

    # WARNING:
    # We may need to resort atoms in the structure corrector function
    # In such case, bonds and charges must be resorted as well and saved apart to keep values coherent
    # Bonds are calculated during the structure corrector but atom charges must be extracted no
    charges = get_charges(filtered_topology_file)
    self.project.get_charges._set_parent_output(self.project, charges)

    print(' * Correcting structure')

    # Set output filenames for the already filtered structure and trajectory
    corrected_structure_filepath = f'{output_directory}/{CORRECTED_STRUCTURE}'
    corrected_structure_file = File(corrected_structure_filepath)
    corrected_trajectory_filepath = f'{output_directory}/{CORRECTED_TRAJECTORY}'
    corrected_trajectory_file = File(corrected_trajectory_filepath)

    # Correct the structure
    # This function reads and or modifies the following MD variables:
    #   snapshots, reference_bonds, register, cache, mercy, trust
    structure_corrector(
        structure = provisional_structure,
        input_trajectory_file = imaged_trajectory_file,
        input_topology_file = filtered_topology_file,
        output_structure_file = corrected_structure_file,
        output_trajectory_file = corrected_trajectory_file,
        MD = self,
        pbc_selection = provisional_pbc_selection,
        snapshots = snapshots,
        register = self.register,
        mercy = self.project.mercy,
        trust = self.project.trust,
        guess_bonds = self.project.guess_bonds
    )

    # If the corrected output exists then use it
    # Otherwise use the previous step files
    # Corrected files are generated only when changes are made in these files
    corrected_structure_file = corrected_structure_file if corrected_structure_file.exists else imaged_structure_file
    corrected_trajectory_file = corrected_trajectory_file if corrected_trajectory_file.exists else imaged_trajectory_file

    # Set for every type of file (structure, trajectory and topology) the input, the last processed step and the output files
    input_and_output_files = [
        (input_structure_file, corrected_structure_file, output_structure_file),
        (input_trajectory_files[0], corrected_trajectory_file, output_trajectory_file),
        (input_topology_file, filtered_topology_file, output_topology_file)
    ]
    # Set a list of intermediate files
    intermediate_files = set([
        converted_structure_file, converted_trajectory_file,
        filtered_structure_file, filtered_trajectory_file,
        imaged_structure_file, imaged_trajectory_file,
    ])
    # Now we must rename files to match the output file
    # Processed files remain with some intermediate filename
    for input_file, processed_file, output_file in input_and_output_files:
        # If the processed file is already the output file then there is nothing to do here
        # This means it was already the input file and no changes were made
        if processed_file == output_file: continue
        # There is a chance that the input files have not been modified
        # This means the input format has already the output format and it is not to be imaged, fitted or corrected
        # However we need the output files to exist and we dont want to rename the original ones to conserve them
        # In order to not duplicate data, we will setup a symbolic link to the input files with the output filepaths
        if processed_file == input_file:
            # If output file exists and its the same as the input file, we can not create a symlink from a file to the same file
            if output_file.exists:
                output_file.remove()
            output_file.set_symlink_to(input_file)
        # Otherwise rename the last intermediate file as the output file
        else:
            # In case the processed file is a symlink we must make sure the symlink is not made to a intermediate step
            # Intermediate steps will be removed further and thus the symlink would break
            # If the symlinks points to the input file there is no problem though
            if processed_file.is_symlink():
                target_file = processed_file.get_symlink()
                if target_file in intermediate_files:
                    target_file.rename_to(output_file)
                else:
                    processed_file.rename_to(output_file)
            # If the files is not a symlink then simply rename it
            else:
                processed_file.rename_to(output_file)

    # Save the internal variables
    self._structure_file = output_structure_file
    self._trajectory_file = output_trajectory_file
    self.project._topology_file = output_topology_file

    # If the input and output file have the same name then overwrite input cksums
    # Thus we avoid this task to run forever
    # DANI: Esto es provisional, hay que prohibir que los inputs se llamen como los outputs
    if input_structure_file == output_structure_file:
        task.cache_cksums['input_structure_file'] = get_cksum_id(output_structure_file)
    if len(input_trajectory_files) == 1 and input_trajectory_files[0] == output_trajectory_file:
        task.cache_cksums['input_trajectory_files'] = get_cksum_id([ output_trajectory_file ])
    if input_topology_file == output_topology_file:
        task.cache_cksums['input_topology_file'] = get_cksum_id(output_topology_file)

    # --- Definitive PBC selection ---

    # Now that we have the corrected structure we can set the definitive PBC atoms
    # Make sure the selection is identical to the provisional selection
    if self.pbc_selection != provisional_pbc_selection:
        raise InputError('PBC selection is not consistent after correcting the structure. '
            'Please consider using a different PBC selection. '
            'Avoid relying in atom distances or elements to avoid this problem.')

    # --- RUNNING FINAL TESTS ------------------------------------------------------------

    # Note that some tests have been run already
    # e.g. stable bonds is run in the structure corrector function

    # Note that tests here do not modify any file

    # Check the trajectory has not sudden jumps
    self.is_trajectory_integral()

    # Make a final test summary
    self.print_tests_summary()

    # Issue some warnings if failed or never run tests are skipped
    self._issue_required_test_warnings()
        
    # --- Cleanup intermediate files

    # Set a list of input files to NOT be removed
    inputs_files = set([ input_structure_file, *input_trajectory_files, input_topology_file ])
    # We must make sure an intermediate file is not actually an input file before deleting it
    removable_files = intermediate_files - inputs_files
    # Now delete every removable file
    for removable_file in removable_files:
        # Note that a broken symlink does not 'exists'
        if removable_file.exists or removable_file.is_symlink():
            removable_file.remove()

    # We must leave something in the directory
    # Otherwise the task will find out it is empty and will delete it
    # And we don't want to delete it or it will be done again in the next run
    # To solve this leave some symbolic file
    if all(self.register.tests.values()):
        manifest_filepath = f'{output_directory}/manifest.txt'
        with open(manifest_filepath, 'w') as file:
            file.write('Input files were processed sucessfully.\n' +
                'Intermediate files were removed to save space.\n'
                'Please do not remove this file\n')