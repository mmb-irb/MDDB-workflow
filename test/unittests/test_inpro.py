import shutil
import pathlib
from mddb_workflow.mwf import Project, MD
from mddb_workflow.utils.constants import TRAJECTORY_INTEGRITY_FLAG


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent/'data'
dummy_dir = data_dir/'input/dummy'
test_fld = data_dir/'output/test_inpro'

def regenerate_test_fld():
    # Remove old test and copy files
    shutil.rmtree(test_fld, ignore_errors=True)
    test_fld.mkdir(parents=True, exist_ok=True)

def test_amber_prmtop_nc():
    regenerate_test_fld()

    # Copy necessary files
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"amber/ala_ala.prmtop", test_fld)
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"amber/trajectory.nc", test_fld/"replica_1")

    # Initialize Project and process files
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.prmtop"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/trajectory.nc"),
                    md_directories=['replica_1'])
    md = project.mds[0]
    md.input_files_processing(md)

def test_amber_top_nc():
    regenerate_test_fld()
    # Copy necessary files
    shutil.copy(dummy_dir/"inputs.yaml", test_fld)
    shutil.copy(dummy_dir/"amber/ala_ala.prmtop", test_fld/"ala_ala.top")
    (test_fld/"replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir/"amber/trajectory.nc", test_fld/"replica_1")

    # Initialize Project and process files
    project = Project(directory=str(test_fld),
                    input_topology_filepath=str(test_fld/"ala_ala.top"),
                    input_trajectory_filepaths=str(test_fld/"replica_1/trajectory.nc"),
                    md_directories=['replica_1'])
    md = project.mds[0]
    md.input_files_processing(md)

# def test_mercy():
#     regenerate_test_fld()
#     shutil.copy(dummy_dir/"inputs.yaml", test_fld)
#     shutil.copy(dummy_dir/"gromacs/ala_ala.tpr", test_fld)
#     (test_fld/"replica_1").mkdir(exist_ok=True)
#     shutil.copy(dummy_dir/"gromacs/trajectory.xtc", test_fld/"replica_1")

#     for mercy in [True, False]:
#         print(f"\nTesting mercy={mercy}\n")
#         # Initialize Project and process files
#         project = Project(directory=str(test_fld),
#                         input_topology_filepath=str(test_fld/"ala_ala.tpr"),
#                         #input_trajectory_filepaths=str(test_fld/"replica_1/trajectory.xtc"),
#                         md_config=[['replica_1','*xtc']],
#                         mercy=mercy)
#         md: MD = project.mds[0]
#         md.register.add_warning(TRAJECTORY_INTEGRITY_FLAG, 'Dummy error to trigger mercy')
#         md.register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
#         md.input_files_processing(md)
    

if __name__ == "__main__":
    test_amber_prmtop_nc()
    test_amber_top_nc()
    # test_mercy()