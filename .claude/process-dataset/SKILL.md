---
name: process-dataset
description: >-
  Process a collection of MDDB projects with the MDDB workflow using the library's Dataset
  class and its `mwf dataset` CLI. Use this whenever the user wants to register many project/MD
  directories into a dataset, generate per-project inputs.yaml files from a template, launch the
  workflow across a dataset (locally, in parallel, or on SLURM), check/monitor the state of
  projects (new/waiting/running/error/done), read their logs, or upload finished projects to an
  MDDB node with `mwf ds load`. Trigger this for any request phrased as "run the workflow over my
  dataset", "add these projects to the dataset", "generate inputs for all projects", "change this
  field in every inputs file", "which projects failed", "re-run the errored ones", "upload the done
  projects", "batch process these MDs", or anything involving the mwf_ds.sql / *.db dataset file,
  even if the user doesn't name the Dataset class explicitly.
---

# Processing a dataset with the MDDB workflow

The `Dataset` class (`mddb_workflow/core/dataset.py`) manages a **collection** of MDDB projects and
their MDs (replicas/subprojects) so the workflow can be run over all of them and their state tracked
centrally. State lives in a single SQLite file (conventionally `mwf_ds.sql` or `*.db`); each project
and MD is a row with a `state` (`new` → `waiting` → `running` → `error`/`done`), a `message`, and a
relative path. UUIDs stored in per-directory `.mwf_cache.json` files let the dataset track entries
even if directories move.

Two equivalent interfaces exist, and choosing correctly matters:

- **The `mwf dataset` CLI** — the default. Use it for interactive, one-off, and operational tasks
  (registering directories, generating inputs, launching, checking status). It is the interface the
  users know and what most requests map to.
- **The `Dataset` Python API** — use it only when the task needs logic the CLI can't express: custom
  filtering, post-processing the status DataFrame, scripting a multi-step pipeline, or writing an
  inputs generator. Don't reach for Python when a single CLI command does the job.

`dataset` has the alias `ds`, so `mwf ds ...` works everywhere `mwf dataset ...` does. The
`dataset_path` argument is **positional and optional** on every subcommand — omit it and the CLI uses
the first `*.db` or `*.sql` file in the current directory. Prefer passing it explicitly when more than
one could exist, so the user always knows which dataset they're touching. Only `add` and `scan` create
a missing dataset file; the other subcommands fail instead.

The CLI changes often. Before relying on a flag, check it with `mwf ds <subcommand> -h` in the
environment that will run it, and note which workflow version (`git log -1` in the workflow repo) a
dataset was processed with. Different machines (a laptop, the cluster) may have different versions.

## The end-to-end flow

Most dataset work is some subset of these six steps, in order. Figure out which step the user is on
rather than always running the whole thing.

1. **Register** projects and MDs into the dataset (`add` or `scan`).
2. **Generate** a per-project `inputs.yaml` from a Jinja2 template (`inputs`), or edit one field
   across them (`inputs -s`).
3. **Launch** the workflow over the selected projects (`run`).
4. **Monitor** progress and inspect state (`show`, `ls`, `logs`, `watch`, `status`, `view`).
5. **Re-run** the ones that failed, after fixing the cause.
6. **Upload** the finished projects to an MDDB node (`load`).

### 1. Register projects and MDs

Two ways to populate the dataset:

- **`add`** — you know the directory patterns. Give glob patterns for projects and, with `-md`, the
  subdirectory patterns that are MDs (replicas) inside each project.

  ```bash
  mwf dataset add mwf_ds.sql -p 'project_*' 'special_cases/*' -i '*/logs' -md 'replica*'
  ```

  `-p/--paths_or_globs` project dirs/globs · `-i/--ignore_dirs` dirs to skip · `-md/--md_dirs` MD
  subdir patterns matched *within* each project.

  **`-md` is usually optional — leave it out by default.** The MDs get registered automatically the
  first time the workflow runs in each project (the run reads that project's `inputs.yaml`, discovers
  its MDs from the `mds` list, and adds their rows). So a plain `mwf dataset add -p '...'` normally
  suffices. Pass `-md` only when you specifically want the per-replica rows in the dataset *before* any
  run — e.g. to inspect or track MD-level status up front.

- **`scan`** — the directories already have `.mwf_cache.json` files (e.g. an existing dataset being
  re-indexed, or dirs that moved). Walks the tree, reads UUIDs from the cache, and adds projects
  (no `project_uuid` in cache) and MDs (has `project_uuid`) accordingly.

  ```bash
  mwf dataset scan .
  ```

  Prefer `scan` when re-registering known projects — it preserves existing UUIDs and paths. Use `add`
  for first-time registration from patterns.

### 2. Generate inputs.yaml per project

Each project needs an `inputs.yaml` describing its metadata and MDs before the workflow runs.
Generate them from one Jinja2 **template** plus an optional **generator** that supplies per-project
values.

```bash
mwf dataset inputs mwf_ds.sql -it template.yaml -ig generator.py -o -p 'project_*' -st new
```

- `-it/--inputs_template` — path to the Jinja2 template. The canonical full template lives at
  `mddb_workflow/resources/inputs_file_template.yml`; copy and trim it for a given dataset.
- `-ig/--inputs_generator` — path to a `.py` file defining `inputs_generator(project_dir)` that
  returns a dict of extra template variables (see `references/inputs-generation.md`). Optional.
- `-o/--overwrite` — regenerate even where `inputs.yaml` already exists. Without it, existing files
  are left untouched — safe by default, so pass `-o` only when the user wants to replace them.
- The query filters — restrict which projects get an inputs file (see **Filtering** below).

The same subcommand reads and edits existing inputs files, which is how to change one setting across a
dataset without regenerating everything:

```bash
mwf dataset inputs mwf_ds.sql -g pbc_selection                      # print one field per project
mwf dataset inputs mwf_ds.sql -s mercy '[stabonds]' -st error -ms '*stable bonds*'
mwf dataset inputs mwf_ds.sql -p 'project_1'                       # without -it: print whole files
```

`-s FIELD VALUE` parses the value as YAML and rewrites the files, which drops their comments. Run the
same filters with `-g` first so the user sees what will change.

The template always receives `DIR` (project dir name), `PATH` (full path), `DATASET` (relative path
to the dataset file), and `project_status` (the project's DB row), plus whatever the generator
returns. Read `references/inputs-generation.md` before writing a template or generator.

Before generating inputs for a large dataset, **dry-run the generator over every project** (fast,
read-only, no `mwf`) to catch half-copied dirs, missing trajectories, or unexpected replica counts
before they become failed jobs — see the "Pre-flight" section in `references/inputs-generation.md`.

**Always set `dataset_path: {{DATASET}}` in the template.** This is the single most important line:
it wires each generated `inputs.yaml` back to its dataset file (as a path relative to the project), so
the workflow updates the right dataset's state when it runs. A template that omits it — or hardcodes a
path — will produce projects that don't report status correctly. `assets/inputs_template.yaml` is a
complete, working example (a GROMACS nanobody dataset) you can copy and adapt; it sets `dataset_path`
and pairs with `assets/nanobodies_generator.py`.

### 3. Launch the workflow

```bash
# Local, sequential
mwf dataset run mwf_ds.sql -p 'project_*' -st new -c 'mwf run -filt -fit'

# Local, 4 projects in parallel
mwf dataset run mwf_ds.sql -st new -ps 4 -c 'mwf run'

# SLURM (requires a job template, not a -c command)
mwf dataset run mwf_ds.sql -st new -sl -jt job_template.sh
```

- `-c/--cmd` — the command run in each project directory (typically an `mwf run ...` invocation).
- `-sl/--slurm` — submit via `sbatch` instead of running locally. **SLURM mode requires `-jt` and
  rejects `-c`** — the job template must itself contain the workflow command.
- `-jt/--job-template` — Jinja2-templated bash/SLURM script rendered per project (receives the same
  `DIR`/`PATH`/`project_status` plus all `inputs.yaml` fields and a `slurm` flag).
  `assets/job_template.sh` is a ready-to-use example. See `references/launching.md`. With neither
  `-c` nor `-jt`, `run` reuses the `mwf_slurm_job.sh` a previous run left in each project, and launches
  nothing if any selected project lacks one.
- `-ps/--pool_size` — number of concurrent local jobs (>1 enables parallel mode; ignored under SLURM).
- `-n/--n_jobs` — cap on how many projects to launch (0 = all). Handy for a test run on a few.
- `--debug` — print what *would* run (and render job scripts) without executing or submitting.
  **Always offer a `--debug` dry run first** for any large or SLURM launch — it's the cheapest way to
  catch a bad filter or template before spawning hundreds of jobs.

`run` operates on **projects only** (not MD rows). The command runs with the project directory as its
working directory.

**Scaling to many jobs — SLURM array.** For very large datasets, submitting one `sbatch` per project
via `run -sl` is heavy on the scheduler. The alternative is a single **SLURM array job** that indexes
into the list of project directories itself, with a concurrency cap (`--array=0-N%K`). This is the
"send a lot of jobs together" pattern. `assets/job_array.sh` is a working example (a nanobody dataset,
220 dirs, 30 concurrent) — submit it directly with `sbatch job_array.sh`, not through `mwf dataset run`. See
`references/launching.md` for how to adapt the array size, directory glob, and the per-task command.

**Workflow flags for large datasets.** Before launching thousands of projects, time one project end to
end with the intended flags, then propose two flags to the user:

- **`-e energies clusters pockets`**: these analyses dominate the runtime (on PLAS-20K `clusters` took
  ~23 s per replica and `energies` minutes, against ~1 min for the rest of the run) and `pockets` writes
  large `.dx` files. Excluding them means they won't be in the database, so confirm with the user and
  record it in the README.
- **`-lb` (local BLAST)**: remote NCBI BLAST drops connections under load and tasks hang until their
  time limit. `-lb` runs `blastp` against a local Swiss-Prot copy. Check the database exists before
  submitting an array, or every task will start downloading it at once.

Details, the database check and the remaining network limits (HTTP 429 from PubChem/RCSB) are in
"Workflow flags for large datasets" in `references/launching.md`.

**On the IRB cluster**, read `references/irb-cluster.md` before writing any `sbatch` script: conda
activation without `module`, using the lab's own partitions (e.g. `mmb_cpu_zen3`/`mmb_cpu_sphr`),
`spot_cpu` (preemptible, requeue, cannot be combined with lab partitions in one job), and the
6000-task array limit.

### 4. Monitor and inspect

```bash
mwf dataset show mwf_ds.sql -m                    # count per state + the 10 most frequent errors
mwf dataset show mwf_ds.sql                       # table of the last 50 entries
mwf dataset show mwf_ds.sql -st done -n 0         # all done rows (0 = no row limit)
mwf dataset ls mwf_ds.sql -st error -n 5          # one block per project: full message, name,
                                                  #   whether it is loaded, MDs
mwf dataset logs mwf_ds.sql -e -st error -n 3     # print the latest error logs
mwf dataset ls mwf_ds.sql -wn -n 0                # projects with any current warning
mwf dataset watch mwf_ds.sql                       # live-updating view
mwf dataset status mwf_ds.sql -p project_1         # one project's status
mwf dataset view mwf_ds.sql -p project_1           # open its processed structure and trajectory in VMD
```

`show` truncates long paths and messages, so use `ls` to read them. `show` flags: `-m/--summary` ·
`-n/--n_rows` row cap (default 50, 0 = all) · `-l/--include_logs` adds links to latest `.out`/`.err` ·
`-s/--sort_by` column · `-r/--reverse` · `-id/--show_id` reveals UUIDs. All the query filters below
apply too. Flags of every subcommand are in `references/cli-reference.md`.

`show`, `ls`, `watch` and `status` auto-reconcile stale `running` rows: if a tracked PID or SLURM job is
gone, the entry is flipped to `error` with an explanatory message. So a project stuck at `running`
after its job died will correct itself the next time you view the dataset. This means these
"read" commands write to the dataset; see the guardrails before inspecting a dataset you must not
change.

**Warnings from the checks.** The workflow stores the warnings of its checks (e.g. `pbcontacts`
for contacts across the box, `intrajrity` for the trajectory integrity) in the project cache, and
`-wn <tags>` selects projects by them. After a run, list the projects with warnings and look at a
few with `view` before uploading: imaging problems are much cheaper to catch here than in the web
client.

### 5. Re-run failures

The state model makes iteration cheap: inspect failures, fix the cause, then launch again filtered to
the errored (or specific) projects.

```bash
mwf dataset show mwf_ds.sql -st error -ms 'URLError*'   # group/inspect a failure class
mwf dataset run mwf_ds.sql -st error -c 'mwf run'        # re-launch only the failed ones
```

**Re-processing a whole dataset with a newer workflow.** To redo every step (e.g. after a fix in the
imaging), pass `-ow` alone to `mwf run`: it re-runs every task and overwrites its outputs. Before
that, check that the raw inputs named in each `inputs.yaml` still exist (scratch cleanups delete them;
see `references/irb-cluster.md`), since `-ow` starts again from them. Projects already uploaded keep
their database entry through the `.project_id` file in the project directory (see step 6). `-ow` only
replaces workflow outputs, so that file survives; don't delete or recreate project directories to
start over.

### 6. Upload finished projects

```bash
mwf dataset ls mwf_ds.sql -st done -n 0                         # check what would be uploaded
mwf dataset load mwf_ds.sql -ld ~/repo/MDDB/loader -nd <node>   # first time: loader dir and node
mwf dataset load mwf_ds.sql -p 'project_1'                      # later: both are remembered
mwf dataset load mwf_ds.sql -o -st done                         # overwrite data already uploaded
```

⚠️ **Outward-facing**: it writes to a remote database, so confirm the selection and get the user's
go-ahead first. `load` selects `done` projects by default. The loader writes a `.project_id` file into
each project it uploads and reads it on the next load, so re-uploading updates the same database
entry instead of creating a new one. **Keep those files**: if one is lost, the next load creates a
duplicate project. `ls` shows which projects have one (`Loaded: yes`). Details, the nodes config and
partial uploads are in `references/data-lifecycle.md`.

### Keep a README booking file

Every dataset gets a `README.md` at its root (next to the `*.db`). It records how the dataset was
processed, so the next person, or the next session, can continue without re-deriving anything.
Datasets often come back months later as a "continuation", and the hard-won facts are otherwise lost
in scripts and old shells. Create it at the start of the work, and update it whenever a decision is
made or a problem is found. Don't wait until the end. If a previous batch left a README (check the
backup dir too), read it first and continue it.

Record the facts that can't be recovered by reading the scripts:

- **Source:** where the raw data lives (remote/path, access command), its layout, and which files
  are actually needed versus skipped, with sizes.
- **Selection:** how many projects were selected and why. Include what was excluded (overlap with
  earlier batches, missing files) and the files listing the selected IDs.
- **Layout and scripts:** the directory layout and what each script in `utils/` does.
- **Differences from earlier batches:** for example new file formats or a changed ligand selection.
- **Commands:** the exact `mwf run` flags and the `mwf dataset` commands per step.
- **Known issues:** per-project problems such as truncated replicas, missing topologies or special
  mercy flags.
- **Log:** a dated list of what was done (download started or finished, array job IDs, error
  classes and their fixes).

## Filtering (query flags)

`inputs`, `show`, `ls`, `logs`, `watch`, `run` and `load` all accept the same filters, and combining
them precisely is what makes dataset operations safe — you target exactly the projects you mean:

- `-p/--query_path` — glob(s) on relative path, e.g. `-p 'project_[1-3]' '*special*'`. Default `*`.
- `-st/--query_state` — one or more states: `new waiting running error done`.
- `-ms/--query_message` — glob(s) on the status message, e.g. `-ms 'URLError*' '*permission*'` —
  useful for isolating one class of failure.
- `-qi/--query_inputs FIELD VALUE` — a field of the inputs file equals a value (YAML-parsed; for a list
  field, any item). Repeat for more conditions, e.g. `-qi temp 310 -qi mds.0.name 'replica 1'`.
- `-nm/--query_name` — glob(s) on the project `name` in its inputs file, ignoring case.
- `-wn/--query_warnings` — glob(s) on the tags of current warnings (MDs included); alone, any warning.
- `-sc/--query_scope` — `project`/`p` or `md`/`m` to restrict scope (`show` and `watch` only).

Filters combine with AND across flag types and OR within a flag's list. When in doubt about what a
filter will hit, run `ls` or `show` with the same filters first — it's the same selection `run` and
`load` will use.

## Choosing the CLI vs. the Python API

Reach for Python (`from mddb_workflow.core.dataset import Dataset, State`) only when the CLI can't
express the task. Typical triggers: post-processing `ds.get_dataframe(...)` with pandas, driving a
multi-step pipeline in one script, or conditional logic across steps. The API mirrors the CLI:
`add_entries` / `scan` → `generate_inputs_yaml` / `set_inputs_field` → `launch` → `get_dataframe` /
`summary` / `error_summary` / `list_projects` / `show_logs` → `load`. Full method signatures and a
worked example are in `references/python-api.md`.

## Guardrails

- **Confirm which dataset file** you're operating on before mutating (`add`, `remove`, `inputs -o`,
  `inputs -s`, `run`, `load`). Passing `dataset_path` explicitly avoids silently hitting the wrong
  `*.db`.
- **"Read-only" commands write.** Opening a dataset updates its schema, and `show`/`ls`/`watch`/
  `status` flip stale `running` rows to `error`. When the user asks for an assessment without changes,
  or the dataset belongs to a run with another workflow version, copy the file and query the copy
  read-only with Python's `sqlite3` (`file:<copy>?mode=ro`).
- **Check the database file itself on shared filesystems.** Many concurrent array tasks writing one
  SQLite file on BeeGFS can leave `OperationalError: disk I/O error` messages and a corrupted index.
  On a copy, run `PRAGMA integrity_check`. If a `sqlite_autoindex_*` index is missing rows (a symptom:
  `count(*)` disagrees with `GROUP BY state`, or a project stays `running` forever), `REINDEX` on the
  real file fixes it. Ask before running it.
- **`load` is outward-facing.** Confirm the selection with `ls` and get a go-ahead. Keep the
  `.project_id` files: losing one makes the next load create a duplicate project.
- **Dry-run big or SLURM launches** with `--debug` and show the user the count and sample commands
  before the real run.
- **Pre-flight the generator** over the whole dataset before generating inputs at scale, and surface
  any projects that fail to resolve or fall in an unexpected replica-count bucket.
- **Never widen a filter to make something "work."** If `run -st error` selects nothing, the fix is to
  inspect state with `show`, not to drop the filter and re-run everything.
- **`remove` deletes rows, not directories** — but it's still destructive to the dataset's record and
  cascades to a project's MDs. Confirm the `-p` pattern with `show` first.
- **Update the dataset's `README.md`** after each step: decisions, counts, job IDs and new problems.
  See "Keep a README booking file".

## Reference files

Read these when the task goes beyond the summary above:

- `references/cli-reference.md` — every `mwf dataset` subcommand with all flags and examples.
- `references/inputs-generation.md` — writing the Jinja2 template and the `inputs_generator` function.
- `references/launching.md` — local/parallel/SLURM execution, job-template authoring, and workflow
  flags for large datasets (`-e` slow analyses, `-lb` local BLAST, rate limits).
- `references/irb-cluster.md` — IRB cluster specifics: environment activation, mmb and spot partitions
  (priority, QoS, preemption), array size limit and chaining, submission pitfalls.
- `references/python-api.md` — the `Dataset` class API and a full programmatic walkthrough.
- `references/data-lifecycle.md` — the steps around processing: downloading raw data (including
  S3/MinIO buckets with rclone: list once, select files, download on the cluster), uploading
  finished projects with `mwf ds load` (the `.project_id` file, nodes config, partial uploads), and
  archiving them (the last two driven by the `done` state).

## Bundled files

**`assets/`** — templates/inputs fed into `mwf`/`sbatch` (copy and adapt per dataset; examples are from
a GROMACS nanobody dataset):

- `inputs_template.yaml` — Jinja2 inputs template following `inputs_schema.py`; sets
  `dataset_path: {{DATASET}}` and builds the `mds` list. Copy and edit the metadata for a new dataset.
- `nanobodies_generator.py` — the matching `inputs_generator(project_dir)`: discovers the `r*/md.xtc`
  replica trajectories, the representative `r1/md.tpr` topology, and the `*_solv_ions.gro` structure,
  returning `mds` + `top` + `struct` for the template.
- `job_template.sh` — per-project SLURM job template for `mwf dataset run -jt` (runs
  `mwf run -img -fit -e energies clusters pockets -lb` on the mmb partitions).
- `job_array.sh` — standalone SLURM **array** script for submitting many projects in one job, with a
  `done`-state guard so re-submitting skips finished projects (same flags and partitions).

**`scripts/`** — standalone executable helpers for the steps *around* processing, run directly (see
`references/data-lifecycle.md`):

- `download_data.py` — acquire raw project archives via `aria2c` and unpack them into project dirs.
- `download_rclone.sh` — SLURM job that downloads a `--files-from` list from an S3/MinIO bucket with
  rclone, keeps files safe from the scratch cleanup, skips excluded projects and lists missing files.
- `upload_projects.py` — legacy: publishes `done` projects via the loader. Use `mwf ds load` instead;
  keep the script only for its continuous `--watch` mode (⚠️ outward-facing).
- `zip_projects.py` — archive `done` projects and reclaim disk (⚠️ deletes originals after zipping).
