"""Tests of the provenance read from the BioBB RO-Crates."""
import json

import pytest

from mddb_workflow.tools.provenance import find_crate_files, get_provenance
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.file import File
from mddb_workflow.utils.quickhash import QUICKHASH_SAMPLE_SIZE, quickhash

pytestmark = [pytest.mark.CI, pytest.mark.unit_int]

RECORDED = 'Skipped by the restart property: the outputs come from a previous execution that was not recorded.'


class Crate:
    """Build a small RO-Crate with the same shape as the BioBB ones."""
    def __init__(self, workflow: str = None):
        self.graph = [{'@id': './', '@type': 'Dataset'}]
        if workflow:
            self.graph.append({'@id': '#workflow-run', '@type': 'CreateAction', 'name': f'Run of {workflow}'})

    def file(self, name: str, size: int, value: str) -> dict:
        file_id = f'file:///{name}-{value}'
        if not any(entity['@id'] == file_id for entity in self.graph):
            self.graph.append({'@id': f'#quickhash-{value}', '@type': 'PropertyValue', 'name': 'quickhash',
                'algo': 'xxh64-sample', 'version': 1, 'sample_size': QUICKHASH_SAMPLE_SIZE, 'value': value})
            self.graph.append({'@id': file_id, '@type': 'File', 'name': name, 'contentSize': str(size),
                'quickhash': {'@id': f'#quickhash-{value}'}})
        return {'@id': file_id}

    def step(self, step_id: str, tool: str, inputs: list, outputs: list, parameters: dict = {},
             recorded: bool = True, mdp: dict = None, agent: str = None, environment: dict = None):
        objects = list(inputs)
        for name, value in parameters.items():
            self.graph.append({'@id': f'#{step_id}-{name}', '@type': 'PropertyValue', 'name': name, 'value': value})
            objects.append({'@id': f'#{step_id}-{name}'})
        if mdp:
            references = []
            for key, value in mdp.items():
                self.graph.append({'@id': f'#{step_id}-mdp-{key}', '@type': 'PropertyValue', 'name': key, 'value': value})
                references.append({'@id': f'#{step_id}-mdp-{key}'})
            self.graph.append({'@id': f'#{step_id}-grompp.mdp', '@type': 'PropertyValue', 'name': 'grompp.mdp',
                'value': f'{len(mdp)} parameters', 'valueReference': references})
            objects.append({'@id': f'#{step_id}-grompp.mdp'})
        if not any(entity['@id'] == f'#{tool}' for entity in self.graph):
            self.graph.append({'@id': f'#{tool}', '@type': 'SoftwareApplication', 'name': tool, 'softwareVersion': '5.3.2'})
        action = {'@id': f'#{step_id}', '@type': 'CreateAction', 'name': f'Run of {tool}',
            'actionStatus': 'http://schema.org/CompletedActionStatus', 'instrument': {'@id': f'#{tool}'},
            'object': objects, 'result': outputs[0] if len(outputs) == 1 else outputs,
            'startTime': '2026-11-19T14:23:21+01:00', 'endTime': '2026-11-19T14:23:56+01:00'}
        if agent:
            self.graph.append({'@id': f'#{agent}', '@type': 'Person', 'name': agent})
            action['agent'] = {'@id': f'#{agent}'}
        if environment:
            references = []
            for key, value in environment.items():
                self.graph.append({'@id': f'#{step_id}-env-{key}', '@type': 'PropertyValue', 'name': key, 'value': value})
                references.append({'@id': f'#{step_id}-env-{key}'})
            self.graph.append({'@id': f'#{step_id}-env', '@type': 'PropertyValue', 'name': 'execution environment',
                'value': 'a summary', 'valueReference': references})
            action['executionEnvironment'] = {'@id': f'#{step_id}-env'}
        else:
            # Old BioBB versions write just a link
            action['executionEnvironment'] = 'https://biobb-common.readthedocs.io/en/latest/provenance.html#executionEnvironment'
        if recorded:
            action['description'] = RECORDED
            del action['startTime']
        self.graph.append(action)

    def write(self, directory):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / 'ro-crate-metadata.json').write_text(json.dumps({'@graph': self.graph}))


def real_file(crate: Crate, path) -> dict:
    """Add a file which exists in disk, so the MD files can be found by their cksum."""
    return crate.file(path.name, path.stat().st_size, quickhash(str(path)))


@pytest.fixture
def project(tmp_path):
    """Two MDs sharing the topology and the structure, plus redundant and branch steps."""
    md_files = {}
    for name, text in [('topology.tpr', 'top'), ('structure.pdb', 'pdb'), ('r1.xtc', 'traj 1'), ('r2.xtc', 'traj 2')]:
        md_files[name] = tmp_path / name
        md_files[name].write_text(text)
    provenance = tmp_path / 'provenance'
    # Setup workflow, shared by both MDs
    setup = Crate('Setup')
    raw = setup.file('raw.pdb', 10, 'aaaa')
    system = setup.file('system.gro', 20, 'bbbb')
    imaged = setup.file('system_img.gro', 20, 'bbbb')
    setup.step('build', 'biobb Build', [raw], [system], {'binary_path': '/home/someone/gmx', 'force_field': 'amber99sb-ildn'})
    setup.step('image', 'biobb GMXImage', [system], [imaged])
    setup.step('tpr', 'biobb Grompp', [system], [real_file(setup, md_files['topology.tpr'])],
        mdp={'dt': 0.002, 'ref-t': '300 300'})
    setup.step('pdb', 'biobb ToPdb', [system], [real_file(setup, md_files['structure.pdb'])])
    setup.step('check', 'biobb Check', [system], [setup.file('check.txt', 5, 'cccc')])
    setup.write(provenance / 'setup')
    # The crate of a single step repeats the step of the workflow crate
    single = Crate()
    single.step('build', 'biobb Build', [single.file('raw.pdb', 10, 'aaaa')], [single.file('system.gro', 20, 'bbbb')])
    single.write(provenance / 'setup' / 'provenance' / 'Build_1')
    # One production workflow per MD
    for number in [1, 2]:
        production = Crate(f'Production {number}')
        topology = real_file(production, md_files['topology.tpr'])
        trajectory = real_file(production, md_files[f'r{number}.xtc'])
        # The first MD declares who ran it, the second only has the user in the machine
        production.step(f'run{number}', 'biobb Mdrun', [topology], [trajectory], recorded=number == 2,
            agent='Jane Doe' if number == 1 else None,
            environment={'username': 'jdoe', 'hostname': 'node1.cluster', 'cpu_count': 64, 'slurm_job_id': '123'})
        production.step(f'view{number}', 'biobb TrjConv', [trajectory], [production.file('view.xtc', 30, f'dddd{number}')])
        production.write(provenance / f'rep{number}')
    md_files = {f'rep{number}': {
        'topology': File(str(md_files['topology.tpr'])),
        'structure': File(str(md_files['structure.pdb'])),
        'trajectories': [File(str(md_files[f'r{number}.xtc']))],
    } for number in [1, 2]}
    return find_crate_files(str(provenance)), md_files


def by_id(part: dict) -> dict:
    return {step['id']: step for step in part['steps']}


def test_steps_are_split_between_project_and_mds(project):
    crate_files, md_files = project
    provenance = get_provenance(crate_files, md_files)
    project_steps = by_id(provenance['project'])
    assert set(project_steps) == {'build', 'image', 'tpr', 'pdb', 'check'}
    assert set(by_id(provenance['mds']['rep1'])) == {'run1', 'view1'}
    assert set(by_id(provenance['mds']['rep2'])) == {'run2', 'view2'}
    # Files of the common steps are only in the project part
    project_files = {file['id'] for file in provenance['project']['files']}
    md_files_ids = {file['id'] for file in provenance['mds']['rep1']['files']}
    assert md_files_ids.isdisjoint(project_files)
    assert provenance['mds']['rep1']['md_files']['topology'] in project_files


def test_hidden_steps_are_marked(project):
    crate_files, md_files = project
    provenance = get_provenance(crate_files, md_files)
    steps = {**by_id(provenance['project']), **by_id(provenance['mds']['rep1'])}
    # All outputs are identical to some input
    assert steps['image'].get('redundant') is True
    # They do not lead to any MD file
    assert steps['check'].get('branch') is True
    assert steps['view1'].get('branch') is True
    for step_id in ['build', 'tpr', 'pdb', 'run1']:
        assert 'redundant' not in steps[step_id] and 'branch' not in steps[step_id]


def test_step_fields(project):
    crate_files, md_files = project
    provenance = get_provenance(crate_files, md_files)
    steps = by_id(provenance['project'])
    build = steps['build']
    # Read from the workflow crate, not the crate of the single step
    assert build['workflow'] == 'Setup'
    assert build['parameters'] == {'force_field': 'amber99sb-ildn'}
    assert build['status'] == 'recorded'
    assert 'start' not in build and 'end' not in build
    assert steps['tpr']['parameters']['grompp.mdp'] == {'dt': 0.002, 'ref-t': '300 300'}
    run1 = by_id(provenance['mds']['rep1'])['run1']
    assert run1['status'] == 'completed'
    assert run1['start'] and run1['end']
    # Files with the same content are the same file, with all their names
    system = next(file for file in provenance['project']['files'] if file['id'] == '20-bbbb')
    assert system['names'] == ['system.gro', 'system_img.gro']


def test_user_and_environment(project):
    crate_files, md_files = project
    provenance = get_provenance(crate_files, md_files)
    run1 = by_id(provenance['mds']['rep1'])['run1']
    run2 = by_id(provenance['mds']['rep2'])['run2']
    # The declared person is preferred over the username in the machine
    assert run1['user'] == 'Jane Doe'
    assert run2['user'] == 'jdoe'
    # The username is not repeated in the environment and the job id is not kept
    assert run1['environment'] == {'hostname': 'node1.cluster', 'cpu_count': 64}
    # Steps with just the link to the documentation have no environment
    build = by_id(provenance['project'])['build']
    assert 'environment' not in build and 'user' not in build


def test_md_file_missing_in_provenance(project, tmp_path):
    crate_files, md_files = project
    other = tmp_path / 'other.xtc'
    other.write_text('not in the provenance')
    md_files['rep2']['trajectories'] = [File(str(other))]
    with pytest.raises(InputError) as excinfo:
        get_provenance(crate_files, md_files)
    assert 'other.xtc' in str(excinfo.value)


def test_no_provenance(tmp_path):
    assert find_crate_files(str(tmp_path / 'provenance')) == []
    assert get_provenance([], {}) is None
