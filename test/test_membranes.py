import os, pytest
from mddb_workflow.mwf import Project
from mddb_workflow.utils.file import File
from mddb_workflow.utils.type_hints import *
from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step


def get_reduced_trajectory(project: Project, test_proj_dir, n_frames: int = 100):
    mdframes = project.remote.get_project_data()["metadata"]["mdFrames"]
    traj_path = f'{test_proj_dir}/replica_1/trajectory.xtc'
    # if os.path.exists(traj_path):
    #     os.remove(traj_path)
    os.makedirs(os.path.dirname(traj_path), exist_ok=True)
    output_file = File(traj_path)
    frame_step, reduced_frame_count = calculate_frame_step(mdframes, n_frames)
    print(os.getcwd(), output_file.path)
    project.remote.download_trajectory(output_file, 
                                frame_selection=f'1:-1:{frame_step}',
                                format='xtc')

class TestMembranes:
    """Tests for membrane-related analyses:
    'membs': ['memmap', 'density',  'thickness', 'apl', 'lorder', 'linter']"""

    @pytest.fixture(scope="class", params=["A01IP", "A01MA", "A01M9"])
    def test_accession(self, request):
        return request.param
    
    def test_thickness(self, project: 'Project', test_proj_dir: str):
        """Test thickness analysis"""
        get_reduced_trajectory(project, test_proj_dir, n_frames=100)
        project.mds[0].overwritables.add('thickness')
        project.mds[0].run_thickness_analysis(project.mds[0])

    def test_lorder(self, project: 'Project', test_proj_dir: str):
        """Test lipid order analysis"""
        get_reduced_trajectory(project, test_proj_dir, n_frames=100)
        project.mds[0].overwritables.add('lorder')
        project.mds[0].run_lipid_order_analysis(project.mds[0])