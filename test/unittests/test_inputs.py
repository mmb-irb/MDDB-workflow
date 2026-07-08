import pytest

from mddb_workflow.core.inputs_schema import (
    validate_inputs,
    KNOWN_INPUT_FIELDS,
)
from mddb_workflow.utils.auxiliar import InputError


@pytest.mark.unit_int
def test_minimal_inputs_are_valid():
    """A minimal inputs file with just a name should validate."""
    data = {'name': 'My project'}
    # Should not raise and should return the very same dict
    assert validate_inputs(data) is data


@pytest.mark.unit_int
def test_empty_inputs_are_valid():
    """An inputs file with no fields set should still validate (all optional)."""
    assert validate_inputs({}) == {}


@pytest.mark.unit_int
def test_full_inputs_are_valid():
    """A rich inputs file exercising the nested structures should validate."""
    data = {
        'name': 'PD-1 in complex with nivolumab',
        'authors': ['Alice', 'Bob'],
        'groups': 'Some group',
        'framestep': 0.01,
        'timestep': 2,
        'temp': 310,
        'type': 'trajectory',
        'ff': ['CHARMM27'],
        'pdb_ids': ['5GGR'],
        'forced_references': {'A': 'noref', 'B': 'Q15116'},
        'chainnames': {'A': 'Antibody', 'C': 'Protein'},
        'ligands': [
            {'pubchem': 1986},
            {'drugbank': 'DB00945'},
            {'chembl': 'CHEMBL25', 'vmd_selection': ['chain D']},
        ],
        'interactions': [
            {
                'name': 'Protein-antibody interaction',
                'agent_1': 'Protein',
                'selection_1': 'chain C',
                'agent_2': 'Antibody',
                'selection_2': 'chain A or chain B',
                'distance_cutoff': 10,
            }
        ],
        'links': [{'name': 'Data source', 'url': 'https://example.org', 'ping': True}],
        'customs': [
            {
                'name': 'A custom view',
                'representations': [
                    {'name': 'lig', 'selection': 'VIR', 'type': 'ball+stick', 'color': 'element'},
                ],
            }
        ],
        'mds': [
            {'name': 'replica 1', 'mdir': 'replica_1', 'temp': 310},
            {'name': 'replica 2', 'mdir': 'replica_2', 'input_trajectory_filepaths': 'traj.xtc'},
        ],
        'mdref': 0,
    }
    assert validate_inputs(data) is data


@pytest.mark.unit_int
def test_forced_references_accepts_list_and_dict():
    """Forced references may be given either as a list or as a per-chain dict."""
    assert validate_inputs({'forced_references': ['Q9BYF1', 'P0DTC2']})
    assert validate_inputs({'forced_references': {'A': 'Q9BYF1', 'C': 'noref'}})


@pytest.mark.unit_int
def test_ligand_requires_an_accession():
    """A ligand without any database accession must be rejected."""
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'ligands': [{'vmd_selection': 'chain D'}]})
    assert 'at least one' in str(excinfo.value)


@pytest.mark.unit_int
def test_ligand_shorthand_pubchem_cid():
    """A bare PubChem CID (int or numeric string) is expanded into a full entry."""
    data = {'ligands': [1986, '2244', {'drugbank': 'DB00945'}]}
    validate_inputs(data)
    assert data['ligands'][0] == {'pubchem': '1986'}
    assert data['ligands'][1] == {'pubchem': '2244'}
    assert data['ligands'][2] == {'drugbank': 'DB00945'}


@pytest.mark.unit_int
def test_ligand_name_string_is_rejected():
    """A non-numeric bare string (a ligand name) must be rejected."""
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'ligands': ['aspirin']})
    assert 'aspirin' in str(excinfo.value)


@pytest.mark.unit_int
def test_interaction_missing_field_is_rejected():
    """An interaction missing a required agent field must be rejected."""
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'interactions': [{'name': 'x', 'agent_1': 'a', 'selection_1': 's'}]})
    message = str(excinfo.value)
    # The error message should point at the offending location
    assert 'interactions' in message
    assert 'agent_2' in message or 'selection_2' in message


@pytest.mark.unit_int
def test_md_type_must_be_valid():
    """The MD 'type' must be either 'trajectory' or 'ensemble'."""
    assert validate_inputs({'type': 'trajectory'})
    assert validate_inputs({'type': 'ensemble'})
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'type': 'md'})
    assert 'trajectory' in str(excinfo.value)


@pytest.mark.unit_int
def test_numeric_fields_must_be_positive():
    """framestep, timestep and temp must be positive when provided."""
    assert validate_inputs({'framestep': 0.01, 'timestep': 2, 'temp': 310})
    for field in ('framestep', 'timestep', 'temp'):
        with pytest.raises(InputError) as excinfo:
            validate_inputs({field: -1})
        assert field in str(excinfo.value)
        with pytest.raises(InputError):
            validate_inputs({field: 0})


@pytest.mark.unit_int
def test_pdb_ids_format_is_validated():
    """PDB ids must match the expected format, as a string or a list."""
    assert validate_inputs({'pdb_ids': '6M0J'})
    assert validate_inputs({'pdb_ids': ['5GGR', '6ACS']})
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'pdb_ids': ['not-a-pdb']})
    assert 'not-a-pdb' in str(excinfo.value)


@pytest.mark.unit_int
def test_wrong_type_is_rejected():
    """A scalar where a list of structured items is expected must be rejected."""
    with pytest.raises(InputError):
        validate_inputs({'interactions': 'not a list'})


@pytest.mark.unit_int
def test_bad_framestep_is_rejected():
    """A non-numeric framestep must be rejected."""
    with pytest.raises(InputError):
        validate_inputs({'framestep': 'not a number'})


@pytest.mark.unit_int
def test_orientation_must_have_16_numbers():
    """The viewer orientation must be a list of exactly 16 numbers."""
    assert validate_inputs({'orientation': list(range(16))})
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'orientation': [1, 2, 3]})
    assert '16' in str(excinfo.value)


@pytest.mark.unit_int
def test_unknown_field_warns_but_passes(capsys):
    """Unknown top-level fields are preserved and only warned about."""
    data = {'name': 'p', 'membranes': None, 'totally_made_up': 42}
    assert validate_inputs(data) is data
    captured = capsys.readouterr()
    assert 'Unknown field' in captured.out
    assert 'membranes' in captured.out
    assert 'totally_made_up' in captured.out

@pytest.mark.unit_int
def test_unknown_field_strict_raises():
    """In strict mode (e.g. forced CLI inputs) unknown fields must raise."""
    with pytest.raises(InputError) as excinfo:
        validate_inputs({'name': 'p', 'totally_made_up': 42}, strict_unknown=True)
    message = str(excinfo.value)
    assert 'totally_made_up' in message
    # The error should hint at the available inputs
    assert 'Available inputs' in message
    assert 'framestep' in message
    # Known fields must still pass in strict mode
    assert validate_inputs({'name': 'p', 'temp': 310}, strict_unknown=True)


@pytest.mark.unit_int
def test_error_message_lists_all_problems():
    """The formatted error should report every problem found, not just the first."""
    data = {
        'framestep': 'not a number',
        'ligands': [{'vmd_selection': 'x'}],
    }
    with pytest.raises(InputError) as excinfo:
        validate_inputs(data)
    message = str(excinfo.value)
    assert 'framestep' in message
    assert 'ligands' in message


@pytest.mark.CI
@pytest.mark.unit_int
def test_known_fields_match_template():
    """Every field declared in the schema must appear in the inputs template
    (and vice versa) so the template and the validator stay in sync."""
    import yaml
    from mddb_workflow.utils.constants import INPUTS_TEMPLATE_FILEPATH

    with open(INPUTS_TEMPLATE_FILEPATH) as handle:
        template = yaml.safe_load(handle)
    template_fields = set(template.keys())
    # Fields the schema knows but the template does not document
    missing_from_template = KNOWN_INPUT_FIELDS - template_fields
    assert not missing_from_template, (
        f'Schema fields missing from the template: {sorted(missing_from_template)}')
    # Fields the template documents but the schema does not know about
    missing_from_schema = template_fields - KNOWN_INPUT_FIELDS
    assert not missing_from_schema, (
        f'Template fields missing from the schema: {sorted(missing_from_schema)}')


@pytest.mark.unit_int
def test_real_example_inputs_validate():
    """The example inputs files shipped with the repo must validate."""
    import glob
    import os
    import yaml

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    patterns = [
        os.path.join(repo_root, 'test', 'data', 'input', '**', 'inputs.yaml'),
        os.path.join(repo_root, 'docs', 'source', 'new_dataset', '**', 'inputs.yaml'),
    ]
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern, recursive=True))
    assert files, 'No example inputs files were found to validate'
    for filepath in files:
        with open(filepath) as handle:
            data = yaml.safe_load(handle)
        # Empty files (e.g. the dummy project) have nothing to validate
        if not data:
            continue
        # Should not raise
        validate_inputs(data)
