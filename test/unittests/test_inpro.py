import os
import shutil
import pathlib
from mddb_workflow.mwf import Project, MD, Structure


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent/'data'
dummy_dir = data_dir/'input/dummy'
test_fld = data_dir/'output/test_inpro'


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
                    md_directories=['replica_1'])
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
                    md_directories=['replica_1'])
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
    Structure.check_incoherent_bonds = lambda _: True
    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)
    for mercy in [True, False]:
        print(f"\nTesting mercy={mercy}\n")
        # Initialize Project and process files
        project = Project(
            input_topology_filepath=str("ala_ala.tpr"),
            md_config=[['replica_1', '*xtc']],
            mercy=mercy)
        md: MD = project.mds[0]
        md.input_files_processing(md)
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
                    md_directories=['replica_1'])
    project.membrane_map
    # Second run to test cache loading
    print("\n ------ Second run to test cache loading ------ \n")
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.tpr"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/raw_trajectory.xtc"),
                    md_directories=['replica_1'])
    project.membrane_map
    os.chdir(cwd)


if __name__ == "__main__":
    test_amber_prmtop_nc()
    test_amber_top_nc()
    # test_mercy()
    # test_mda_universe()
