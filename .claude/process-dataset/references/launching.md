# Launching the workflow over a dataset

`mwf dataset run` (or `Dataset.launch`) runs a command for **each project** matched by the query
filters. Three execution modes.

## 1. Local sequential (default)

Runs one project at a time in its own directory, streaming output to both the terminal and per-run
log files under `<project>/logs/mwf_<timestamp>.out` / `.err`.

```bash
mwf dataset run mwf_ds.sql -st new -c 'mwf run -filt -fit'
```

## 2. Local parallel

Set `-ps/--pool_size` > 1 to run several projects concurrently (capped at the CPU count). Output goes
to the log files (not streamed live to the terminal), with progress printed as jobs complete.

```bash
mwf dataset run mwf_ds.sql -st new -ps 4 -c 'mwf run'
```

## 3. SLURM

Submits each project via `sbatch`. **SLURM mode requires `-jt/--job-template` and rejects `-c`** — the
workflow command must live inside the job template, because the template is what SLURM executes.
Submitted projects are marked `waiting` until the job starts.

```bash
mwf dataset run mwf_ds.sql -st new -sl -jt job_template.sh
```

## The job template (`assets/job_template.sh`)

A bash/SLURM script rendered with Jinja2 once per project. It receives:

- `DIR` — project directory name.
- `PATH` — full project path.
- `project_status` — the project's DB row dict.
- All key/values from that project's `inputs.yaml` (so you can branch on project metadata).
- `slurm` — boolean, true when launched with `-sl`.

The rendered script is written to `<project>/mwf_slurm_job.sh` and made executable. Example:

```bash
#!/bin/bash
#SBATCH --partition=mmb_cpu_zen3,mmb_cpu_sphr
#SBATCH --job-name=mwf_{{DIR}}
#SBATCH --output=logs/mwf_%j.out
#SBATCH --error=logs/mwf_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/miniforge3/etc/profile.d/conda.sh   # `module` is not available on the IRB compute nodes
conda activate mwf_env
mwf run -filt -fit -e energies clusters pockets -lb
```

On the IRB cluster, partitions, QoS, spot preemption and array limits are in `irb-cluster.md`.

**Reusing the previous job script.** If `run` gets neither `-c` nor `-jt`, it reuses the
`mwf_slurm_job.sh` rendered in each project by an earlier `-jt` launch. All selected projects must
have one, otherwise nothing is launched and the missing ones are listed. Handy to re-submit failures
with exactly the same settings:

```bash
mwf dataset run mwf_ds.sql -st error -ms 'Workflow process is no longer running*' -sl
```

SLURM submission uses `--output=logs/mwf_%j.out --error=logs/mwf_%j.err` with the project directory as
the working directory, so the `logs/` directory is created automatically.

## Scaling: the SLURM array job (`assets/job_array.sh`)

For datasets with hundreds or thousands of projects, one `sbatch` per project (`run -sl`) floods the
scheduler. A single **SLURM array job** is the better tool: it enumerates the project directories once
and lets SLURM run them as array tasks with a built-in concurrency cap. This is what "send a lot of
jobs together" means. Submit it directly — **not** through `mwf dataset run`:

```bash
sbatch job_array.sh
# or override the range/concurrency at submit time:
sbatch --array=0-50%2 job_array.sh
```

`assets/job_array.sh` is a working example. The parts to adapt for a new dataset:

- **`#SBATCH --array=0-N%K`** — `N` = (number of project dirs − 1), `K` = max concurrent tasks. Get the
  count with e.g. `ls -d vhh_run/*/ | wc -l` (the example is `0-219%30` for 220 projects).
- **The directory list** — `DIRS=($(ls -d vhh_run/*/))` then `DIR=${DIRS[$SLURM_ARRAY_TASK_ID]}`. Point
  the glob at your project directories, filtering out any non-project entries (e.g. append
  `| grep -vE '(\.zip$|\.aria2$)'` if the folder also holds archives). The array range must match this
  list's length.
- **The per-task command** — `cd "$DIR"` then the `mwf run ...` invocation (here
  `mwf run -img -fit -e energies clusters pockets -lb`). Set the flags for your analyses; see "Workflow
  flags for large datasets" below.
- **Optional per-task steps** (the state guard is active in the example; the inputs step is commented):
  - `mwf dataset status ../../nanobodies.db -p .` to read the current state and `exit 0` early if
    already `done` — makes re-submitting the array cheap since finished projects skip immediately.
  - `mwf dataset inputs ../../nanobodies.db -it $UTILS/inputs_template.yaml -ig $UTILS/nanobodies_generator.py
    -p "*${DIR_NAME}" -o` to (re)generate that project's `inputs.yaml` just before running (skip if you
    generated all inputs in bulk beforehand).
- **Log handling** — array logs land in `array_logs/mwf_%A_%a.out/.err` at the submit dir and are then
  copied into each project's `logs/`. Create `array_logs/` before submitting.

Trade-off vs. `run -sl`: the array job is far lighter on the scheduler and easy to resume, but it
bypasses the dataset's automatic `waiting`/`running` bookkeeping that `mwf dataset run` provides — so
lean on the `mwf dataset status` guard inside the script (and `mwf dataset show` afterward) to track
state.

## Workflow flags for large datasets

Flags that are harmless on ten projects decide whether thousands finish. Settle them on one project
before launching anything wide.

**Time one project end to end first.** Run a single representative project interactively (on IRB with
`srun`, see `irb-cluster.md`) with the intended flags, then rank the tasks by time:

```bash
sed 's/\x1b\[[0-9;]*m//g' run.log | grep "completed in" | awk '{print $2, $5}' | sort -k2 -rn | head
```

This also surfaces blocking input problems in minutes instead of after a whole array has failed.

**Exclude the slow analyses with `-e`.** On PLAS-20K (6.4k atoms, 5 × 40 frames per project)
`clusters` alone took ~23 s per replica and `energies` ran for minutes before being cancelled, while a
full run without them took ~1 min. `pockets` also wrote ~5 MB `.dx` files per pocket (~60 MB per
replica). Unless the user needs them, exclude them for datasets in the thousands:

```bash
mwf run -fit -e energies clusters pockets ...
```

Excluding is a user decision (the analyses will be missing from the database), so propose it with the
timings and confirm. Record the excluded list in the dataset README. When re-running a project that
already ran with these analyses, delete their leftovers (`rep*/energies`, `rep*/incomplete_energies`,
`rep*/clusters`, `rep*/pockets`, `rep*/incomplete_pockets`): an interrupted energies run left a
192 MB `restart` file and made the per-project disk estimate ten times too high.

**Use local BLAST with `-lb`.** The protein mapping runs a BLAST for every chain without a UniProt
reference. The remote NCBI service often drops connections ("Remote end closed connection without
response"), and after the 3 retries the task can hang until its time limit. With many concurrent tasks
this happens every run. `-lb/--local_blast` runs `blastp` locally against Swiss-Prot instead; on
PLAS-20K it turned projects that hung for 20+ min into 1–5 min runs. Requirements:

- BLAST+ (`blastp`, `update_blastdb.pl`) in the environment (`conda install -c bioconda blast`).
- The Swiss-Prot database in `mddb_workflow/resources/blastdb/`. It is downloaded on first use, so
  **check it exists before submitting an array**, otherwise every task starts downloading it at once:
  ```bash
  python -c "from mddb_workflow.utils.blast import local_blastdb_exists, BLASTDB_DIRECTORY; \
             print(BLASTDB_DIRECTORY, local_blastdb_exists())"
  ```
  If it prints `False`, run one project with `-lb` first (or `update_blastdb.pl --decompress swissprot`
  in that directory) and launch the rest afterwards.
- The local copy is a snapshot (check the file dates), while remote NCBI is live, so a few UniProt
  matches may differ. Note the database date in the README.

`-lb` only replaces BLAST. Other external requests in the workflow (PubChem, the RCSB ligand page,
InterProScan) still go over the network and return HTTP 429 when too many tasks hit them at once from
the same IP. Keep the whole run's concurrency moderate and add a random start delay to each array task
(`sleep $((RANDOM % 61))` before `mwf run`), so that tasks do not all reach the ligand step together.

## Workflow flags for re-processing and imaging

These `mwf run` flags (or the matching inputs file fields) matter when re-processing a dataset with a
newer workflow, typically to fix trajectories that look wrong in the web client. Check
`mwf run -h` for the installed version.

- **`-ow` alone** re-runs every task and overwrites its outputs. `-ow <tasks...>` re-runs only those.
  It starts again from the raw inputs, so check they still exist first.
- **`-img`** images the trajectory; **`-fit [auto]`** fits it (`auto`: only if it does not look fitted
  already); **`-trans auto|X Y Z`** translates the system before imaging (`auto`, the default in the
  inputs template, finds a translation that keeps non-PBC atoms away from the box boundaries). The
  inputs fields `image`, `fit` and `translation` can also be set per MD inside `mds`.
- **`-pbc <selection>` / `pbc_selection`** decides which atoms stay under PBC after imaging.
  `auto` guesses solvent, counter ions **and lipids**. Centring, the molecule-aware no-jump and the
  automatic translation only act on the atoms *outside* this selection. So in a lipid-only system
  (bilayers with water and ions) `auto` leaves nothing to centre on, and imaging is just a
  `-pbc res`. If the bilayer must be centred or kept whole across the box, set the selection to the
  non-lipid part, e.g. `mwf ds inputs mwf_ds.sql -s pbc_selection 'water or ions'`, and check a
  pilot with `mwf ds view` before the full run.
- **`-t/--trust <checks>`** skips checks; **`-m/--mercy <checks>`** lets them fail without stopping the
  run. Checks: `stabonds`, `cohbonds`, `intrajrity` (trajectory integrity: sudden jumps, permanent
  fragment breaks), `pbcontacts` (contacts across the box); mercy also accepts `elements`, `refseq`,
  `interact`, `ligands`, `lipids`, `chains`, `largeaa`. Failed checks under mercy become warnings that
  `mwf ds ls -wn <tag>` finds afterwards.

## Controlling scope and dry-running

- `-n/--n_jobs N` — launch at most `N` projects (0 = all). Use a small `N` for a pilot run.
- `--debug` — render job scripts and print the commands that *would* run, without executing or
  submitting. **Do this first for any large or SLURM launch** and show the user the project count and
  a sample command before the real run.

```bash
# Dry run: see what a SLURM launch would submit
mwf dataset run mwf_ds.sql -st new -sl -jt job_template.sh --debug

# Pilot: just the first 3 matched projects, locally
mwf dataset run mwf_ds.sql -st new -n 3 -c 'mwf run'
```

## State tracking during a run

`run` records process/SLURM job identifiers on each entry. If a `running` project's process or SLURM
job later disappears, the next `show`/`ls`/`watch`/`status` reconciles it to `error` with a message
noting likely external termination (time limit, OOM, cancellation). A `waiting` project whose job left
the queue without starting gets its own message. This means a crashed run won't stay stuck at
`running` — view the dataset and it self-corrects.

If a project stays `running` although no job exists and viewing does not fix it, check the dataset
file itself (`PRAGMA integrity_check` on a copy): with many array tasks writing one SQLite file on a
shared filesystem, a row can drop out of the `rel_path` index, and then the workflow can no longer
update it. `REINDEX` fixes it (ask first). Use `mwf ds logs -e -p <project>` to read why it failed.
