import json

import pytest

from mddb_workflow.core.inputs import ProjectInputs
from mddb_workflow.mwf import Project, input_files
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.file import File


pytestmark = pytest.mark.CI


def _project_inputs_with_mds(tmp_path, md_config):
    inputs_path = tmp_path / 'inputs.json'
    inputs_path.write_text(json.dumps({'mds': md_config}), encoding='utf-8')
    return ProjectInputs(File(str(inputs_path)), directory_name='demo')


def test_project_inputs_loads_and_resolves_values(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'name: "$DIR project"\n'
        'description: from file\n',
        encoding='utf-8',
    )

    inputs = ProjectInputs(
        File(str(inputs_path)),
        directory_name='demo',
        forced_inputs=[['description', 'from CLI']],
    )

    assert inputs.file_inputs['name'] == 'demo project'
    assert inputs.get('description') == 'from CLI'
    assert inputs.get('name') == 'demo project'


def test_project_passes_raw_forced_inputs_to_project_inputs(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'name: from-file\n'
        'mds:\n'
        '  - name: md1\n'
        '    mdir: md1\n',
        encoding='utf-8',
    )

    project = Project(
        directory=str(tmp_path),
        inputs_filepath=str(inputs_path),
        forced_inputs=[['name', 'from-cli']],
    )

    assert project.inputs.forced_inputs == {'name': 'from-cli'}
    assert project.inputs.get('name') == 'from-cli'
    assert project.get_file_input('name') == 'from-cli'
    assert project.inputs.get_inputs_file() == File(str(inputs_path))
    assert input_files['inputs'](project) == project.inputs.inputs_file
    assert ProjectInputs(File(str(inputs_path)), directory_name='demo').get('name') == 'from-cli'


def test_project_inputs_merges_new_md_configurations(tmp_path):
    inputs = _project_inputs_with_mds(tmp_path, [{'mdir': 'existing', 'name': 'old'}])

    merged = inputs.merge_md_config(
        [['existing', 'existing.xtc'], ['new', 'new.xtc']],
        input_topology_filepath='project.tpr',
        input_structure_filepath='project.pdb',
    )

    expected = [
        {
            'mdir': 'existing',
            'name': 'old',
            'input_topology_filepath': 'project.tpr',
            'input_structure_filepath': 'project.pdb',
            'input_trajectory_filepaths': ['existing.xtc'],
        },
        {
            'mdir': 'new',
            'input_topology_filepath': 'project.tpr',
            'input_structure_filepath': 'project.pdb',
            'input_trajectory_filepaths': ['new.xtc'],
        },
    ]
    assert merged == expected

    reloaded = ProjectInputs(inputs.inputs_file, directory_name='demo')
    assert reloaded.file_inputs['mds'] == expected


def test_project_inputs_normalizes_null_and_removed_md_configurations(tmp_path):
    inputs = _project_inputs_with_mds(tmp_path, [
        {'mdir': 'removed', 'removed': True},
        {'mdir': 'active'},
    ])

    merged = inputs.merge_md_config(None)

    assert merged[0] == {'mdir': 'removed', 'removed': True}
    assert merged[1] == {'mdir': 'active'}

    reloaded = ProjectInputs(inputs.inputs_file, directory_name='demo')
    assert reloaded.file_inputs['mds'] == [
        {'mdir': 'removed', 'removed': True},
        {'mdir': 'active'},
    ]


def test_project_inputs_adds_new_md_after_removed_configurations(tmp_path):
    inputs = _project_inputs_with_mds(tmp_path, [{'mdir': 'removed', 'removed': True}])

    merged = inputs.merge_md_config(
        [['new', 'trajectory.xtc']],
    )

    assert merged[1]['mdir'] == 'new'
    assert merged[1]['input_trajectory_filepaths'] == ['trajectory.xtc']


def test_project_inputs_merges_legacy_md_directories_and_warns(tmp_path, capsys):
    inputs = _project_inputs_with_mds(tmp_path, [{'mdir': 'existing'}])

    merged = inputs.merge_md_config(None, ['existing', 'new'])

    assert merged == [{'mdir': 'existing'}, {'mdir': 'new'}]
    assert 'deprecated' in capsys.readouterr().out


def test_project_inputs_rejects_duplicate_legacy_md_directories(tmp_path):
    inputs = _project_inputs_with_mds(tmp_path, [{'mdir': 'existing'}])

    with pytest.raises(InputError, match='duplicated MD directories'):
        inputs.merge_md_config(None, ['new', 'new'])

    reloaded = ProjectInputs(inputs.inputs_file, directory_name='demo')
    assert reloaded.file_inputs['mds'] == [{'mdir': 'existing'}]


def test_project_inputs_merges_and_persists_md_config(tmp_path, capsys):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text('mds: null\n', encoding='utf-8')
    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    resolved = inputs.merge_md_config(None, ['new'])

    assert resolved == [{'mdir': 'new'}]
    assert 'deprecated' in capsys.readouterr().out
    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.file_inputs['mds'] == [{'mdir': 'new'}]


def test_project_inputs_expands_legacy_md_directory_globs(tmp_path, monkeypatch):
    (tmp_path / 'replica_1').mkdir()
    (tmp_path / 'replica_2').mkdir()
    monkeypatch.chdir(tmp_path)
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
    )

    merged = inputs.merge_md_config(None, ['replica_*'])

    assert {config['mdir'] for config in merged} == {'replica_1', 'replica_2'}


def test_project_inputs_applies_common_md_overrides_and_skips_removed(tmp_path):
    inputs = _project_inputs_with_mds(tmp_path, [
        {'mdir': 'removed', 'removed': True},
        {'mdir': 'active'},
    ])

    result = inputs.merge_md_config(
        None,
        input_topology_filepath='project.tpr',
        input_structure_filepath='project.pdb',
        input_trajectory_filepaths=['project.xtc'],
    )

    assert result[1] == {
        'mdir': 'active',
        'input_topology_filepath': 'project.tpr',
        'input_structure_filepath': 'project.pdb',
        'input_trajectory_filepaths': ['project.xtc'],
    }



@pytest.mark.parametrize(
    ('md_config', 'expected_message'),
    [
        ([], 'There must be at least one MD'),
        ([{}], 'no name and no directory'),
        (
            [{'name': 'same', 'mdir': 'one'}, {'name': 'same', 'mdir': 'two'}],
            'Duplicated values in MD inputs',
        ),
        (
            [{'name': 'one', 'mdir': 'same'}, {'name': 'two', 'mdir': 'same'}],
            'Duplicated values in MD inputs',
        ),
    ],
)
def test_project_inputs_validates_md_config(md_config, expected_message, capsys):
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
    )

    with pytest.raises(InputError, match=expected_message):
        inputs.validate_md_config(md_config)

    if 'Duplicated' in expected_message:
        output = capsys.readouterr().out
        assert 'same name' in output or 'same directory' in output


def test_project_uses_project_inputs_to_autofill_md_identifiers(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'mds:\n'
        '  - name: production replica\n',
        encoding='utf-8',
    )

    project = Project(directory=str(tmp_path), inputs_filepath=str(inputs_path))

    assert project.md_config == [{'name': 'production replica', 'mdir': 'production_replica'}]
    reloaded = ProjectInputs(File(str(inputs_path)), directory_name=tmp_path.name)
    assert reloaded.file_inputs['mds'][0]['mdir'] == 'production_replica'


@pytest.mark.parametrize(
    ('arg_config', 'expected_topology', 'expected_structure', 'expected_trajectories'),
    [
        (['md1', 'top.tpr', 'trajectory.xtc'], 'top.tpr', None, ['trajectory.xtc']),
        (['md1', 'structure.pdb', 'trajectory.xtc'], None, 'structure.pdb', ['trajectory.xtc']),
        (['md1', 'trajectory.xtc'], None, None, ['trajectory.xtc']),
    ],
)
def test_project_inputs_detects_md_file_roles(
    arg_config,
    expected_topology,
    expected_structure,
    expected_trajectories,
):
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
    )

    merged = inputs.merge_md_config([arg_config])

    assert merged[0]['input_topology_filepath'] == expected_topology
    assert merged[0]['input_structure_filepath'] == expected_structure
    assert merged[0]['input_trajectory_filepaths'] == expected_trajectories


@pytest.mark.parametrize(
    ('input_md_config', 'expected_message'),
    [
        ([['md1']], 'Wrong MD configuration'),
        ([['md1', 'traj1.xtc'], ['md1', 'traj2.xtc']], 'duplicated MD directories'),
    ],
)
def test_project_inputs_rejects_invalid_new_md_configurations(input_md_config, expected_message):
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
    )

    with pytest.raises(InputError) as excinfo:
        inputs.merge_md_config(input_md_config)

    assert expected_message in str(excinfo.value)


def test_project_inputs_maps_legacy_pdb_ids_before_validation(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text('pdbIds: 5GGR\n', encoding='utf-8')

    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    assert inputs.file_inputs['pdb_ids'] == '5GGR'
    assert 'pdbIds' not in inputs.file_inputs


def test_project_inputs_effective_values_are_combined_without_mutating_file(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'name: demo\n'
        'mds:\n'
        '  - name: md1\n'
        '    mdir: md1\n',
        encoding='utf-8',
    )

    inputs = ProjectInputs(
        File(str(inputs_path)),
        directory_name='demo',
        forced_inputs=[['description', 'from CLI']],
    )

    effective = inputs.effective_inputs
    assert effective['name'] == 'demo'
    assert effective['description'] == 'from CLI'
    assert 'pbc_selection' not in effective

    effective['mds'][0]['name'] = 'changed'
    assert inputs.file_inputs['mds'][0]['name'] == 'md1'


def test_project_inputs_applies_forced_inputs_and_reports_changes(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text('name: from-file\n', encoding='utf-8')
    inputs = ProjectInputs(
        File(str(inputs_path)),
        directory_name='demo',
        forced_inputs=[['name', 'from-cli'], ['description', 'generated']],
    )

    assert inputs.apply_forced_inputs() is None
    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.file_inputs == {'name': 'from-cli', 'description': 'generated'}


def test_project_inputs_uses_defaults_without_a_file():
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
    )

    assert inputs.get('pbc_selection') == 'auto'
    assert inputs.get('unknown') is None


def test_project_inputs_uses_explicit_fallback_without_warning(tmp_path, capsys):
    inputs = ProjectInputs(
        File(str(tmp_path / 'missing.yaml')),
        directory_name='demo',
    )

    assert 'Missing inputs file. Allowed tasks will be very limited.' in capsys.readouterr().out
    assert inputs.get('name', 'fallback') == 'fallback'
    assert capsys.readouterr().out == ''


def test_project_inputs_warns_only_once_for_a_default(tmp_path, capsys):
    inputs = ProjectInputs(
        File(str(tmp_path / 'missing.yaml')),
        directory_name='demo',
    )

    assert inputs.get('pbc_selection') == 'auto'
    assert inputs.get('pbc_selection') == 'auto'

    output = capsys.readouterr().out
    assert output.count('Missing input "pbc_selection"') == 1
    assert output.count('Guessing atoms under Periodic Boundary Conditions') == 1


def test_project_inputs_loads_from_a_lazy_input_loader(tmp_path):
    downloaded_path = tmp_path / 'downloaded.yaml'
    downloaded_path.write_text('name: downloaded\n', encoding='utf-8')
    missing_path = tmp_path / 'remote.yaml'
    loader_calls = []

    def load_inputs_file(output_file):
        loader_calls.append(output_file)
        File(str(downloaded_path)).copy_to(output_file)

    inputs = ProjectInputs(
        File(str(missing_path)),
        directory_name='demo',
        input_loader=load_inputs_file,
    )

    assert inputs.is_file_available
    assert not inputs.file_exists
    assert loader_calls == []
    assert inputs.get('name') == 'downloaded'
    assert loader_calls == [File(str(missing_path))]
    assert inputs.file_exists
    assert inputs.inputs_file == File(str(missing_path))
    assert inputs.file_inputs is inputs.file_inputs


def test_project_inputs_rejects_an_empty_file(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text('', encoding='utf-8')

    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    with pytest.raises(InputError, match='Input file is empty'):
        _ = inputs.file_inputs


def test_project_inputs_preserves_unknown_file_fields_with_a_warning(tmp_path, capsys):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'name: demo\n'
        'future_field: preserved\n',
        encoding='utf-8',
    )

    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    assert inputs.file_inputs['future_field'] == 'preserved'
    output = capsys.readouterr().out
    assert 'Unknown field' in output
    assert 'future_field' in output


def test_project_inputs_rejects_invalid_file(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text('type: invalid\n', encoding='utf-8')

    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    with pytest.raises(InputError):
        _ = inputs.file_inputs


def test_project_inputs_validates_known_forced_values():
    with pytest.raises(InputError, match='type'):
        ProjectInputs(
            File('/tmp/project-inputs-does-not-exist.yaml'),
            directory_name='demo',
            forced_inputs=[['type', 'invalid']],
        )


def test_project_inputs_rejects_unknown_forced_input():
    with pytest.raises(InputError):
        ProjectInputs(
            File('/tmp/project-inputs-does-not-exist.yaml'),
            directory_name='demo',
            forced_inputs=[['not_an_input', 'value']],
        )


@pytest.mark.parametrize(
    ('forced_inputs', 'expected_message'),
    [
        ([[]], 'There is an empty "-fin"'),
        ([['name']], 'missing the new input value'),
        ([['name', 'one', 'two']], 'Too many values'),
    ],
)
def test_project_inputs_rejects_malformed_forced_input_pairs(forced_inputs, expected_message):
    with pytest.raises(InputError) as excinfo:
        ProjectInputs(
            File('/tmp/project-inputs-does-not-exist.yaml'),
            directory_name='demo',
            forced_inputs=forced_inputs,
        )

    assert expected_message in str(excinfo.value)


def test_project_inputs_validates_effective_values():
    inputs = ProjectInputs(
        File('/tmp/project-inputs-does-not-exist.yaml'),
        directory_name='demo',
        forced_inputs=[['type', 'trajectory']],
    )

    assert inputs.validate() == {'type': 'trajectory'}
    with pytest.raises(InputError, match='type'):
        inputs.validate({'type': 'invalid'})


@pytest.mark.parametrize(
    ('extension', 'content'),
    [
        ('yaml', 'name: before\n'),
        ('json', '{"name": "before"}\n'),
    ],
)
def test_project_inputs_updates_and_persists_json_and_yaml(tmp_path, extension, content):
    inputs_path = tmp_path / f'inputs.{extension}'
    inputs_path.write_text(content, encoding='utf-8')
    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    assert inputs.update_file_inputs('name', 'after') is True
    assert inputs.update_file_inputs('name', 'after') is False

    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.get('name') == 'after'


def test_project_inputs_updates_nested_list_values(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'mds:\n'
        '  - name: old\n'
        '    mdir: md1\n',
        encoding='utf-8',
    )
    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    assert inputs.update_file_inputs('mds.0.name', 'new') is True

    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.file_inputs['mds'][0]['name'] == 'new'


def test_md_update_keeps_runtime_and_persisted_configs_in_sync(tmp_path):
    inputs_path = tmp_path / 'inputs.json'
    inputs_path.write_text(json.dumps({
        'mds': [{
            'name': 'md1',
            'mdir': 'md1',
            'input_topology_filepath': 'https://example.org/topology.tpr',
        }],
    }), encoding='utf-8')
    project = Project(directory=str(tmp_path), inputs_filepath=str(inputs_path))
    md = project.mds[0]

    assert project.md_config is project.inputs.file_inputs['mds']
    assert md.update_file_inputs('input_topology_filepath', 'topology.tpr') is True

    assert md.get_file_input('input_topology_filepath') == 'topology.tpr'
    assert project.md_config[0]['input_topology_filepath'] == 'topology.tpr'
    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.file_inputs['mds'][0]['input_topology_filepath'] == 'topology.tpr'


def test_project_inputs_update_is_noop_without_a_file(tmp_path):
    inputs = ProjectInputs(
        File(str(tmp_path / 'missing.yaml')),
        directory_name='demo',
    )

    assert inputs.update_file_inputs('name', 'new') is False


def test_project_inputs_creates_the_file_from_a_template(tmp_path, capsys):
    template_path = tmp_path / 'template.yaml'
    inputs_path = tmp_path / 'inputs.yaml'
    template_path.write_text('name: from-template\n', encoding='utf-8')
    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    assert inputs.ensure_file_from_template(File(str(template_path))) is True
    assert inputs_path.exists()
    assert inputs.get('name') == 'from-template'

    assert inputs.ensure_file_from_template(File(str(template_path))) is False
    output = capsys.readouterr().out
    assert 'has been created from a template' in output
    assert 'Inputs file already exists' in output


def test_project_inputs_autofills_values_and_reports_changes(tmp_path):
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(
        'name: demo\n'
        'pbc_selection: old\n',
        encoding='utf-8',
    )
    inputs = ProjectInputs(File(str(inputs_path)), directory_name='demo')

    updates = inputs.autofill({
        'pbc_selection': 'new',
        'description': 'generated',
    })
    assert updates == {'pbc_selection': True, 'description': True}
    assert inputs.get('pbc_selection') == 'new'
    assert inputs.get('description') == 'generated'

    assert inputs.autofill({'pbc_selection': 'new'}) == {'pbc_selection': False}
    reloaded = ProjectInputs(File(str(inputs_path)), directory_name='demo')
    assert reloaded.get('pbc_selection') == 'new'
    assert reloaded.get('description') == 'generated'
