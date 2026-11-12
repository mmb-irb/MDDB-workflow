# Dataset Python API

Use the API only when the CLI can't express the task — custom filtering, pandas post-processing,
multi-step scripting, or an in-process generator. Otherwise prefer `mwf dataset`.

```python
from mddb_workflow.core.dataset import Dataset, State
```

## Lifecycle

```python
ds = Dataset(dataset_path='mwf_ds.sql')               # opens an existing dataset
ds = Dataset(dataset_path='new.db', create=True)      # create it if missing (only add/scan do this in the CLI)
...
ds.close()                                             # release the lock + connection when done
```

`Dataset(...)` raises `InputError` if the file does not exist and `create` is not set. Opening a
dataset runs the schema setup (missing tables, columns and triggers are added), so even opening it
writes to the file. Transient SQLite errors on open are retried with backoff.

Access is guarded by a directory-based lock (`.lock_<db name>`) that is safe on NFS/BeeGFS, so
concurrent CLI and script use of the same dataset is coordinated.

## Query filters

Every method that selects projects takes the same keyword filters as the CLI:

```python
query_path=['project_*'],              # -p   globs on rel_path
query_state=['error'],                 # -st  states
query_message=['URLError*'],           # -ms  globs on the status message
query_inputs=[('temp', '310')],        # -qi  (field, value) pairs in the inputs file, all must match;
                                       #      values parsed as YAML, dots/integers for nesting
query_name=['*kinase*'],               # -nm  globs on the inputs `name`, case ignored
query_warnings=['pbcontacts'],         # -wn  globs on current warning tags ('*' = any warning)
```

`ds.query_table('projects' | 'mds', query_path, query_state, query_message, query_inputs, query_name,
query_warnings)` returns the raw matching rows (tuples in `ds.project_columns` / `ds.md_columns` order)
without reconciling stale `running` rows.

## Registering entries

```python
ds.add_entries(
    paths_or_globs=['project_*', 'special_cases/*'],
    ignore_dirs=['*/logs'],
    md_dirs=['replica*'],
    verbose=True,
)

ds.scan(root_dir='.', ignore_dirs=['scratch'], verbose=True)   # from .mwf_cache.json files
```

## Generating, reading and editing inputs

```python
ds.generate_inputs_yaml(
    inputs_template_path='template.yaml',
    inputs_generator=my_generator,   # a callable(project_dir)->dict, OR a path string to a .py file
    overwrite=False,
    inputs_filename='inputs.yaml',
    query_path=['project_*'],
    query_state=['new'],
)

ds.get_inputs_field('pbc_selection', query_state=['done'])        # prints "<rel_path>: <value>"
ds.set_inputs_field('mercy', '[stabonds]', query_message=['*stable bonds*'])   # value parsed as YAML
ds.show_inputs_files(query_path=['project_1'])                     # prints whole files
```

`set_inputs_field` loads every selected file before writing any (all or nothing), creates missing dict
fields, requires list items to exist, and drops comments from the files it rewrites. See
`inputs-generation.md` for the template variables and generator contract.

## Launching

```python
ds.launch(
    cmd='mwf run -filt -fit',   # OR job_template for SLURM; with neither, the mwf_slurm_job.sh
                                #    of a previous run in each project is reused
    n_jobs=0,                    # 0 = all matched projects
    pool_size=4,                 # >1 → local parallel
    slurm=False,
    job_template=None,           # cmd must be None when slurm=True
    query_path=['project_*'],
    query_state=['new'],
    debug=False,                 # True = print/render only
)
```

## Inspecting state

```python
df = ds.get_dataframe(
    query_path=['project_*'],
    query_state=['error'],
    query_message=['URLError*'],
    query_scope='project',       # 'project'/'p' or 'md'/'m'
    include_logs=True,
    sort_by='last_modified',
    asc=False,
    show_id=False,
)

ds.summary()          # DataFrame: count of projects per state
ds.error_summary()    # DataFrame: error messages grouped by frequency
ds.get_status('path/to/project')   # single entry's status dict

ds.list_projects(n_projects=0, include_logs=True, query_state=['error'])   # = `mwf ds ls`
ds.show_logs(error=True, n_projects=3, query_state=['error'])              # = `mwf ds logs -e`

# Convenience properties
ds.dataframe          # full joined projects+MDs view
ds.projects_table     # raw project rows
ds.mds_table          # raw MD rows
```

`get_dataframe` (so also `list_projects`, `show`, `watch`) and `get_status` reconcile stale `running`
and `waiting` entries against live PIDs and SLURM jobs before returning, flipping dead ones to
`error`. For a strictly read-only look, open a copy of the file with `sqlite3`
(`sqlite3.connect('file:copy.sql?mode=ro', uri=True)`) instead of `Dataset`.

## Uploading

```python
ds.load(
    loader_path='~/repo/MDDB/loader',   # saved in the dataset after the first successful use
    node='<node alias>',                 # from loader_nodes.yaml; also saved
    overwrite=False,                     # True = loader -o, False = loader -c (conserve)
    include=None, exclude=None,          # file globs relative to each project, not both
    skip_trajectories=False, skip_files=False, skip_analyses=False,
    md_directories=None,
    query_state=['done'],                # default
)
```

⚠️ Outward-facing: see `data-lifecycle.md` for the `.project_id` file and when to ask first. It raises
`RuntimeError` at the first project that fails.

## Worked example: re-run everything that hit a specific error

```python
from mddb_workflow.core.dataset import Dataset

ds = Dataset('mwf_ds.sql')
try:
    failed = ds.get_dataframe(query_state=['error'], query_message=['URLError*'])
    print(f"Re-running {len(failed)} projects that failed with a URL error")
    ds.launch(cmd='mwf run', query_state=['error'], query_message=['URLError*'])
finally:
    ds.close()
```

## State enum

```python
State.NEW      # 'new'      — registered, never run
State.WAITING  # 'waiting'  — submitted to SLURM, not yet started
State.RUNNING  # 'running'  — executing (PID / SLURM job tracked)
State.ERROR    # 'error'    — failed (message stored)
State.DONE     # 'done'     — completed successfully
```

## Directory & cache conventions

- Per-directory UUID cache: `.mwf_cache.json` (constant `CACHE_FILENAME` in
  `mddb_workflow/utils/constants.py`). Projects have no `project_uuid` in their cache; MDs do. If a
  cache is lost, the UUID already registered for that path is recovered instead of replacing the
  dataset entry.
- Dataset file: conventionally `mwf_ds.sql` or any `*.db`. Paths in the DB are stored relative to the
  dataset file's directory, so the dataset and its projects can be relocated together.
- Loader trace: `.project_id` in each uploaded project directory (see `data-lifecycle.md`).
