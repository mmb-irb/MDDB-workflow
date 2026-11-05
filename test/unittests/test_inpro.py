import os
import shutil
import pathlib
import pytest
from mddb_workflow.mwf import Project, MD, Structure
from mddb_workflow.utils.auxiliar import TestFailure


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent/'data'
dummy_dir = data_dir/'input/dummy'
test_fld = data_dir/'output/test_inpro'
# The box of the dummy Gromacs trajectory is smaller than the peptide, so it touches its own periodic image
DUMMY_GROMACS_TRUST = ['pbcontacts']


def regenerate_test_fld():
    """Regenerate the test folder by removing the old one and creating a new one."""
    shutil.rmtree(test_fld, ignore_errors=True)
    test_fld.mkdir(parents=True, exist_ok=True)


def test_amber_prmtop_nc():
    """Test processing of Amber prmtop and nc files."""
    regenerate_test_fld()

    # Copy necessary files
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"amber/ala_ala.prmtop", test_fld)
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"amber/raw_trajectory.nc", test_fld/"replica_1")

    # Initialize Project and process files
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.prmtop"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/raw_trajectory.nc"),
                    input_md_directories=['replica_1'])
    md = project.mds[0]
    md.input_files_processing(md)


def test_amber_top_nc():
    """Test processing of Amber top and nc files."""
    regenerate_test_fld()
    # Copy necessary files
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"amber/ala_ala.prmtop", test_fld/"ala_ala.top")
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"amber/raw_trajectory.nc", test_fld/"replica_1")

    # Initialize Project and process files
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.top"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/raw_trajectory.nc"),
                    input_md_directories=['replica_1'])
    md = project.mds[0]
    md.input_files_processing(md)


def test_mercy():
    """Test mercy parameter during trajectory integrity check."""
    regenerate_test_fld()
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"gromacs/ala_ala.tpr", test_fld)
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"gromacs/raw_trajectory.xtc", test_fld/"replica_1")

    # Modify check_incoherent_bonds to force a trajectory integrity failure
    old_fn = Structure.check_incoherent_bonds
    Structure.check_incoherent_bonds = lambda _, ua_selection = None: [0]
    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)
    try:
        # Without mercy the failed test stops the processing
        print("\nTesting mercy=False\n")
        project = Project(
            input_topology_filepath=str("ala_ala.tpr"),
            input_md_config=[['replica_1', '*xtc']],
            mercy=False)
        md: MD = project.mds[0]
        with pytest.raises(TestFailure):
            md.input_files_processing(md)
        # With mercy the processing is completed although the test fails
        # It is run twice so the second time the '*xtc' pattern also finds the output trajectory
        for _ in range(2):
            print("\nTesting mercy=True\n")
            project = Project(
                input_topology_filepath=str("ala_ala.tpr"),
                input_md_config=[['replica_1', '*xtc']],
                mercy=True)
            md: MD = project.mds[0]
            md.input_files_processing(md)
    finally:
        os.chdir(cwd)
        # Restore original function in case there are further tests
        Structure.check_incoherent_bonds = old_fn


def test_mda_universe():
    """Test processing of Gromacs tpr and xtc files using MDAnalysis Universe."""
    regenerate_test_fld()

    # Copy necessary files
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"gromacs/ala_ala.tpr", test_fld)
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"gromacs/raw_trajectory.xtc", test_fld/"replica_1")

    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)
    # Initialize Project and process files
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.tpr"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/raw_trajectory.xtc"),
                    input_md_directories=['replica_1'],
                    trust=DUMMY_GROMACS_TRUST)
    project.membrane_map
    # Second run to test cache loading
    print("\n ------ Second run to test cache loading ------ \n")
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.tpr"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/raw_trajectory.xtc"),
                    input_md_directories=['replica_1'],
                    trust=DUMMY_GROMACS_TRUST)
    project.membrane_map
    os.chdir(cwd)


@pytest.mark.parametrize('ua_selection', ['auto', 'resname ALA'])
def test_ua_selection(ua_selection: str):
    """Test the united atom selection is passed to the coherent bonds test."""
    regenerate_test_fld()
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"gromacs/ala_ala.tpr", test_fld)
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"gromacs/raw_trajectory.xtc", test_fld/"replica_1")

    # Keep the united atom selections passed to check_incoherent_bonds
    passed_selections = []
    old_fn = Structure.check_incoherent_bonds
    def check_incoherent_bonds(structure, ua_selection = None):
        passed_selections.append(ua_selection)
        return old_fn(structure, ua_selection)
    Structure.check_incoherent_bonds = check_incoherent_bonds
    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)
    try:
        project = Project(
            input_topology_filepath=str("ala_ala.tpr"),
            input_md_config=[['replica_1', 'raw_trajectory.xtc']],
            ua_selection=ua_selection,
            trust=DUMMY_GROMACS_TRUST)
        md: MD = project.mds[0]
        md.input_files_processing(md)
        assert len(passed_selections) == 1
        # The dipeptide is atomistic, so the automatic guess finds no united atoms
        if ua_selection == 'auto':
            assert not passed_selections[0]
        else:
            assert passed_selections[0]
            assert passed_selections[0] == md.structure.select(ua_selection)
    finally:
        os.chdir(cwd)
        Structure.check_incoherent_bonds = old_fn


if __name__ == "__main__":
    test_amber_prmtop_nc()
    test_amber_top_nc()
    # test_mercy()
    # test_mda_universe()
