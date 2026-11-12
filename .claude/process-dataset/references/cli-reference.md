# `mwf dataset` CLI reference

Entry point: `mwf` (console script → `mddb_workflow.console:main`). Subcommand group: `dataset`
(alias `ds`). Defined in `mddb_workflow/console.py`; the methods behind it are in
`mddb_workflow/core/dataset.py`. Run `mwf ds <subcommand> -h` to check the flags of the installed
version, since they change often.

`dataset_path` is a **positional, optional** argument on every subcommand. If omitted, the first
`*.db` or `*.sql` file in the current directory is used. Pass it explicitly whenever ambiguity is
possible.

**Only `add` and `scan` create a dataset file.** Every other subcommand raises `InputError` when the
file does not exist, so a typo in the path no longer creates an empty dataset.

**Opening a dataset writes to it.** Every subcommand runs the schema setup when it opens the file
(it adds missing tables, columns and triggers, e.g. the `settings` table used by `load`). `show`,
`watch`, `ls` and `status` also flip stale `running`/`waiting` rows to `error`. To inspect a dataset
strictly read-only (e.g. one another workflow version is still using), copy the file and query the
copy with `sqlite3` in read-only mode (`file:<path>?mode=ro`).

## Subcommands at a glance

| Subcommand | Purpose | Writes? |
|---|---|---|
| `add` | Register projects (and optionally MDs) from glob patterns | DB rows |
| `scan` | Register from existing `.mwf_cache.json` files | DB rows |
| `remove` | Delete project rows (cascades to MDs) | DB rows |
| `status` | One project's or MD's status | stale-row fix |
| `inputs` | Generate, show, read (`-g`) or edit (`-s`) `inputs.yaml` files | inputs files |
| `show` | Table of entries, or a state/error summary (`-m`) | stale-row fix |
| `ls` | One block per project: full path, name, state, message, loaded or not, MDs | stale-row fix |
| `logs` | Print the latest output (or `-e` error) log of each project | no |
| `watch` | Live-updating table | stale-row fix |
| `view` | Open VMD with a project's most processed structure and trajectory | no |
| `run` | Launch the workflow (or any command) per project | DB rows, job scripts |
| `load` | Upload projects to an MDDB node with the loader | ⚠️ remote database |

## `add` — register projects and MDs from patterns

```
mwf dataset add [dataset_path] -p <globs...> [-i <globs...>] [-md <globs...>]
```

| Flag | Long | Meaning |
|------|------|---------|
| `-p` | `--paths_or_globs` | Project directories or glob patterns to register |
| `-i` | `--ignore_dirs` | Directories/globs to skip (default: none) |
| `-md` | `--md_dirs` | Glob patterns for MD (replica) subdirs *within* each project |

```bash
mwf dataset add mwf_ds.sql -p 'project_*' 'special_cases/*' -i '*/logs' -md 'replica*'
```

`-md` is normally omitted: MDs are registered automatically the first time the workflow runs in a
project (it reads the project's `inputs.yaml`, discovers its `mds`, and adds those rows). Use `-md`
only to pre-register per-replica rows before any run. A plain `mwf dataset add -p '...'` is the common
case.

## `scan` — register from existing `.mwf_cache.json` files

```
mwf dataset scan [dataset_path]
```

Walks the current directory tree, reads UUIDs from `.mwf_cache.json`, and adds projects (cache has no
`project_uuid`) and MDs (cache has `project_uuid`). Preserves existing UUIDs and updates paths if a
directory moved. If a cache was lost, the UUID already registered for that path is recovered instead
of replacing the dataset entry.

```bash
mwf dataset scan .
```

## `remove` — delete an entry from the dataset

```
mwf dataset remove [dataset_path] -p <path>
```

`-p/--query_path` selects the project to remove. Removing a project cascades to its MD rows via the
DB foreign key. This deletes **rows**, not files on disk. Confirm the pattern with `show` first.

## `status` — one entry's status

```
mwf dataset status [dataset_path] -p <project_or_md_path>
```

Prints UUID, path, state, scope (`Project`/`MD`), number of MDs, last modified and message. Trailing
slashes and a leading `./` in the path are ignored, so `-p .` works from inside a project.

## `inputs` — generate, show, read or edit inputs files

```
mwf dataset inputs [dataset_path] -it <template> [-ig <generator.py>] [-o] [query filters]   # generate
mwf dataset inputs [dataset_path] [query filters]                                             # show whole files
mwf dataset inputs [dataset_path] -g <field> [query filters]                                  # read one field
mwf dataset inputs [dataset_path] -s <field> <value> [query filters]                          # edit one field
```

| Flag | Long | Meaning |
|------|------|---------|
| `-it` | `--inputs_template` | Jinja2 template path. Without it (and without `-g`/`-s`) the existing inputs files are printed |
| `-ig` | `--inputs_generator` | `.py` file defining `inputs_generator(project_dir)` (needs `-it`) |
| `-o` | `--overwrite` | Regenerate even where inputs.yaml exists (needs `-it`; default: skip existing) |
| `-g` | `--get FIELD` | Print `<rel_path>: <value>` of one field per project. Dots for nesting, integers for list items (`mds.0.name`) |
| `-s` | `--set FIELD VALUE` | Set one field in every selected inputs file. `VALUE` is parsed as YAML (`310` → int, `true`, `[a, b]`, `null`). Missing dict fields are created; list items must exist |

`-g` and `-s` cannot be combined with each other or with `-it/-ig/-o`. `-s` loads every selected file
before writing any, so it either edits all of them or none, but it **drops comments** from the files
it rewrites. Show the selection with `-g <field>` before a `-s`.

```bash
mwf ds inputs nmr.sql -g pbc_selection -p 'BilayerData/*'           # what each project uses now
mwf ds inputs nmr.sql -s pbc_selection 'water or ions' -nm '*POPC*'  # change it for a subset
```

See `inputs-generation.md` for the template variables and the generator contract.

## `show` — display the dataset table

```
mwf dataset show [dataset_path] [-m] [-n N] [-l] [-s COL] [-r] [-id] [query filters] [-sc SCOPE]
```

| Flag | Long | Meaning |
|------|------|---------|
| `-m` | `--summary` | Count of projects per state, then the 10 most frequent error messages (`-n` changes how many) |
| `-n` | `--n_rows` | Max rows (default 50; `0` = all). Shows the last rows after sorting |
| `-l` | `--include_logs` | Add links to latest `.out`/`.err` per entry |
| `-s` | `--sort_by` | Column to sort by (default `last_modified`) |
| `-r` | `--reverse` | Reverse the sort order |
| `-id` | `--show_id` | Reveal `uuid`/`project_uuid` columns |
| `-sc` | `--query_scope` | `project`/`p` or `md`/`m` |

## `ls` — one block per project

```
mwf dataset ls [dataset_path] [-n N] [-l] [-s COL] [-r] [query filters]
```

Prints each project in its own block with the **full** path and message (`show` truncates them), the
`name` from its inputs file, its state, whether it is **loaded** (a `.project_id` file exists in the
project, left there by the loader) and its MDs with their own states and messages. `-l` adds the
paths of the latest logs. `-n` (default 50, `0` = all) keeps the last projects after sorting. Prefer
`ls` over `show` to read error messages or to check which projects are already uploaded.

## `logs` — print the latest logs

```
mwf dataset logs [dataset_path] [-e] [-n N] [query filters]
```

Prints the latest output log (or with `-e/--error` the latest error log) of each selected project, the
most recently modified ones last. `-n` caps the number of projects (default 50, `0` = all). Combine it
with filters to read one failure class: `mwf ds logs -e -st error -ms '*stable bonds*' -n 3`.

## `watch` — live-updating view

```
mwf dataset watch [dataset_path] [-r] [query filters] [-sc SCOPE]
```

## `view` — open a project in VMD

```
mwf dataset view [dataset_path] -p <project_path> [-start S] [-end E] [-step N | -fr F... | -a] [-sel SELECTION]
```

Opens VMD (must be installed where you run it) with the most processed structure and trajectory of
the project, with the PBC box shown. Trajectories under 500 MB load all frames; otherwise pick frames
with `-start/-end/-step`, `-fr` or `-a/--all`. `-sel` (VMD syntax) focuses the view on some atoms and
draws them in CPK. The path must be a project, not an MD directory. Useful to check imaging problems
before uploading, instead of finding them in the web client.

## `run` — launch the workflow over projects

```
mwf dataset run [dataset_path] [-c CMD | -jt TEMPLATE] [-sl] [-ps N] [-n N] [--debug] [query filters]
```

| Flag | Long | Meaning |
|------|------|---------|
| `-c` | `--cmd` | Command to run in each project dir (e.g. `'mwf run -filt -fit'`) |
| `-sl` | `--slurm` | Submit via `sbatch`; rejects `-c` |
| `-jt` | `--job-template` | Jinja2 bash/SLURM script rendered per project. **If neither `-c` nor `-jt` is given, the `mwf_slurm_job.sh` left by a previous run in each project is reused** (all of them must exist, otherwise nothing is launched) |
| `-ps` | `--pool_size` | Concurrent local jobs (>1 → parallel; ignored under SLURM) |
| `-n` | `--n_jobs` | Cap on number of projects launched (0 = all) |
|  | `--debug` | Print the project count and commands/render scripts without executing or submitting |

Operates on **projects only**. See `launching.md`.

## `load` — upload projects to an MDDB node ⚠️

```
mwf dataset load [dataset_path] [-ld LOADER] [-nd NODE] [-o] [-i GLOBS... | -e GLOBS...]
                 [-stj] [-sf] [-sa] [-md DIRS...] [query filters]
```

Outward-facing: it writes to a remote database. Confirm the selection with `ls` (same filters) and
get the user's go-ahead before running it.

| Flag | Long | Meaning |
|------|------|---------|
| `-ld` | `--loader` | Loader directory or its `index.js`. **Saved in the dataset**: pass it once |
| `-nd` | `--node` | Node alias from `loader_nodes.yaml`. **Saved in the dataset**: pass it once |
| `-o` | `--overwrite` | Overwrite data already in the database. Without it only new data is uploaded (loader `-c`) |
| `-i` | `--include_files` | Load only these files (paths/wildcards relative to each project, e.g. `'*/mda.rmsf.json'`) |
| `-e` | `--exclude_files` | Load everything but these. Cannot be combined with `-i` |
| `-stj` | `--skip_trajectories` | Do not load the main trajectory of each MD |
| `-sf` | `--skip_files` | Do not load any file (structures, trajectories...) |
| `-sa` | `--skip_analyses` | Do not load any analysis |
| `-md` | `--md_dirs` | Load only these MD directories (names relative to each project) |

The default state filter of `load` is `-st done`, unlike the other subcommands. The loader's option
counters are updated once at the end instead of after every project. It stops at the first project
that fails. Details (the `.project_id` file, the nodes config, overwriting) are in
`data-lifecycle.md`.

## Query filters (shared by `inputs`, `show`, `ls`, `logs`, `watch`, `run`, `load`)

| Flag | Long | Meaning |
|------|------|---------|
| `-p` | `--query_path` | Glob(s) on relative path (default `*`). Trailing `/` and leading `./` are ignored |
| `-st` | `--query_state` | State(s): `new waiting running error done` |
| `-ms` | `--query_message` | Glob(s) on the status message (e.g. `'URLError*'`) |
| `-qi` | `--query_inputs FIELD VALUE` | A field of the inputs file equals `VALUE` (parsed as YAML; for a list field, one of its items). Dots/integers for nesting (`mds.0.name`). Repeat the flag for more conditions, all must match |
| `-nm` | `--query_name` | Glob(s) on the `name` in the inputs file, ignoring case (e.g. `'*kinase*'`). Projects without a name are left out |
| `-wn` | `--query_warnings` | Glob(s) on the tags of the project's current warnings, MD warnings included (e.g. `pbcontacts`, `intrajrity`). The flag alone selects projects with any warning |
| `-sc` | `--query_scope` | `project`/`p` or `md`/`m` (`show` and `watch` only) |

Filters AND across flag types, OR within a flag's value list.

```bash
mwf ds ls nmr.sql -wn pbcontacts intrajrity -n 0       # projects whose new checks warned
mwf ds show nmr.sql -qi temp 310 -st done              # done projects at 310 K
mwf ds run nmr.sql -nm '*DPPC*' -st error -c 'mwf run -img'
```
