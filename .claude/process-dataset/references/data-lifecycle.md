# Data lifecycle: download, upload, archive

Processing is the middle of a longer pipeline. Uploading is built into the CLI (`mwf ds load`);
downloading and archiving use helper scripts (standalone executables in `scripts/`). Upload and
archive are driven by the dataset's `done` state, which makes the dataset the natural work queue: a
project becomes eligible to publish/archive exactly when the workflow marks it `done`.

```
  download  ─►  add + inputs + run   ─►  mwf ds load    ─►  archive
 (acquire)      (this skill's core)      (publish done)     (zip done, reclaim disk)
```

The scripts are **adaptable examples**, not turnkey tools — each has a CONFIG block at the top with the
dataset-specific paths to set. Upload and archive are outward-facing or irreversible; read the
warnings.

## `download_data.py` — acquire raw projects

Builds a link list from a list of entry IDs, downloads with `aria2c` (resumable, segmented), and
unpacks each archive into a project directory. Requires the `aria2c` binary.

- **Adapt**: the CONFIG block (`API_BASE`, `ID_LIST_PATH`, `DEST_DIR`) and `unzip_project`, which
  encodes *this* source's archive layout (it expects a `protein/` subfolder, falling back to the archive
  root). A different provider means a different `unzip_project`.
- Re-running is safe: entries already downloaded are skipped.
- On scratch, `touch` every downloaded file in the same job (see `irb-cluster.md`): downloads that keep
  their original dates are deleted by the cleanup policy, and then a later `mwf run -ow` has nothing to
  start from.
- After downloading, register what landed with `mwf dataset add` / `scan` and continue with the core
  flow.

## S3/MinIO buckets with rclone

Some datasets arrive as a bucket on an S3-compatible object store (e.g. the JSC MinIO,
`rclone lsd minio_jsc:<bucket>/<path>`). `scripts/download_rclone.sh` is the download job. These are
the lessons from the MDBind batch (2762 projects in the bucket, 979 GB needed out of ~1.7 TB):

**1. Inspect before downloading anything**
- `rclone lsd` gives the project folders. Then `rclone ls <remote>/<one_id>` shows one project's
  layout, and `rclone size` its total.
- `rclone cat <remote>/<id>/<file>.out | sed -n '/&cntrl/,/\//p'` reads the MD parameters (steps,
  dt, write frequency, temperature) without downloading. Check them against the inputs template.
- Look for debris at the bucket root (loose files or stray `preprod_files/`-like folders that
  duplicate project content) and leave it out.

**2. List once, save the listing, derive everything from it**
- A recursive listing is slow: ~20 min for 2764 folders (~370k objects). Do it **once**, in the
  background, and save it:
  ```bash
  rclone lsf -R --files-only --format sp --separator ' ' --checkers 32 \
      --filter '- *.equil.rep*' --filter '+ /*/*.parm7' --filter '+ /*/rep_prod_files/*.rep*.nc' \
      --filter '- *' <remote> > needed_listing.txt
  ```
  Use ordered `--filter` rules; mixing `--include` and `--exclude` warns that their order is
  undefined.
- From that file, compute without touching the bucket again:
  - the total download size;
  - projects missing the topology (MDBind: `1qpb`);
  - the replica count per project;
  - replicas of uneven size within a project, which flag truncated trajectories (MDBind: `1njs`
    rep4 had ~31 of 50 frames).
- Keep only the files the workflow needs (topology + production trajectories). That was ~60% of
  each MDBind project; equilibration, restart and log files were skipped.
- Compare the IDs with earlier batches (old dataset DB, old project lists) and drop the overlap
  before downloading. In MDBind that saved 601 projects and 220 GB.
- Write the result as a `--files-from` list (paths relative to the remote) and keep it in `utils/`.
  It is both the download plan and the record of what was selected.

**3. Download on the cluster, not through a local mount**
- Pushing ~1 TB through the sshfs mount would route it all through the local machine. Run
  `rclone` on the cluster as a SLURM job instead (`irb-cluster.md`: never on the login node).
- The cluster needs the remote. **Ask the user first**, because this copies credentials. Then append
  only that remote's section of the local `rclone.conf` to `~/.config/rclone/rclone.conf` on the
  cluster, `chmod 600`, without printing the keys:
  ```bash
  awk '/^\[/{p=($0=="[<remote>]")} p' "$(rclone config file | tail -1)" | ssh <cluster> \
      'mkdir -p ~/.config/rclone && cat >> ~/.config/rclone/rclone.conf && chmod 600 ~/.config/rclone/rclone.conf'
  ```
- Check that the endpoint is reachable from the cluster: `curl` to it returns HTTP 403 without
  credentials. Then test the remote with `rclone lsf <remote>/<one_id>`. On IRB, compute nodes
  reached the JSC endpoint.
- `rclone copy <remote> <dest> --files-from <list> --no-traverse --transfers 16 --checkers 32` moved
  707 GiB in 35 min (~300 MiB/s). It is resumable, so re-running only fetches what is missing.

**4. Keep the files alive and complete**
- **Use `--local-no-set-modtime --size-only` and `touch` afterwards on scratch.** rclone keeps the
  remote modification time by default, so the scratch cleanup sees fresh downloads as old and deletes
  them: 569 MDBind projects lost their inputs mid-run (see `irb-cluster.md`). The script does all
  three.
- **Make the download repeatable.** Write the missing-files list before and after each run. Filter
  out projects moved to `excluded/`, so a re-download doesn't bring them back and shift the array
  indices.
- **Keep the bucket layout.** Trajectory paths in `inputs.yaml` resolve relative to the project
  directory, so files can stay where they land (e.g. `rep_prod_files/<id>.repN.nc`). There is no
  need to move them into `replica_N/`.
- **Verifying against the bucket costs another listing.** `rclone check --files-from --size-only
  --one-way` lists the bucket again. Prefer the missing-files list plus a per-project count of
  topology and trajectories.
- **Stop remote jobs on the cluster.** Killing a local `ssh` that runs `rclone` doesn't always stop
  it there. Check with `pgrep -af rclone` on the cluster and kill by PID; `pkill -f` can match its own
  shell.

## `mwf ds load` — publish finished projects

```bash
mwf ds ls mwf_ds.sql -st done -n 0                              # what would be uploaded, and
                                                                #   which ones are loaded already
mwf ds load mwf_ds.sql -ld ~/repo/MDDB/loader -nd <node> -p 'project_1'   # pilot on one project
mwf ds load mwf_ds.sql                                          # every done project
mwf ds load mwf_ds.sql -o -p 'project_1'                        # overwrite data already uploaded
mwf ds load mwf_ds.sql -o -i '*/mda.rmsf.json'                  # re-upload only one analysis
```

- ⚠️ **Outward-facing**: it writes to a remote database. Confirm the selection with `ls` (same
  filters), upload one project first, and get the user's go-ahead before loading at scale.
- **Setup.** It needs Node.js, the loader repository (`-ld`, the directory or its `index.js`) and a
  node entry in `loader_nodes.yaml` at the root of the workflow repository: the node alias with
  `ssh_hostname` (an alias in `~/.ssh/config` that forwards the database port, e.g.
  `LocalForward <port> 127.0.0.1:27017`) and the `db_*` fields of that node's loader `.env`. The
  first `load` checks the SSH connection and the database access, then saves `-ld` and `-nd` in the
  dataset (`settings` table), so later loads need neither.
- **Default selection is `done`.** All other query filters apply (`-p`, `-nm`, `-qi`, `-wn`...).
- **Conserve vs overwrite.** Without `-o` the loader runs with `-c`: it uploads what is missing and
  keeps what is there. With `-o` it replaces the data already in the database, which is what a
  re-processed project needs. Narrow it with `-i`/`-e` (file globs relative to each project, not
  combinable), `-md` (only some MD directories), `-stj` (no main trajectories), `-sf` (no files at
  all) or `-sa` (no analyses).
- **It stops at the first failing project.** The option counters are updated once at the end (or
  after the failure), not after each project. Re-running the same command is safe: already-loaded
  data is conserved and the projects keep their database entry.

### The `.project_id` file

When the loader starts on a project, it decides which database entry to use, in this order:

1. an accession forced on the command line or in the inputs file (`accession:`);
2. the id in the project's `.project_id` file, if that project still exists in the database;
3. otherwise, a new project with a new accession.

It then writes the id into `.project_id` in the project directory. So re-uploading a re-processed
project updates the same entry, **as long as `.project_id` is still there**. Consequences:

- **Never delete `.project_id` files**, and don't delete or recreate project directories that hold
  one. If the file is lost, the next load creates a duplicate project. On a scratch filesystem with an
  age-based cleanup, back them up (with their relative paths) outside scratch, or write the accession
  into each `inputs.yaml` (`mwf ds inputs -s accession <A0XXX> -p <project>`) so it no longer depends
  on the file.
- If the id points to a project that no longer exists in the database, the loader warns, removes the
  file and creates a new project.
- The file is written when the load **starts**, so a project with `.project_id` may still be a
  partial upload. Compare with the database (the MDDB tools or the web client) before assuming
  it is complete, and re-run `load` on it.
- `mwf ds ls` shows `Loaded: yes` for projects that have the file.

## `upload_projects.py` — legacy upload script

Superseded by `mwf ds load` above; prefer the CLI. The script is still useful for one thing the CLI
does not do: a continuous worker that uploads projects as they reach `done` while a dataset is still
processing. It polls the dataset for `state=done` projects and uploads each new one via the MDDB
`loader` (`node <loader> load . -c`), recording it in a resume file so re-runs skip it.

- ⚠️ **Outward-facing**: uploading publishes data to an external database. Set `DB_PATH` and the
  `LOADER_CANDIDATES` paths correctly and upload one project first to confirm before running at scale.
- **Two modes.** By default it does one pass (upload everything currently `done`, then exit). Set
  `CONTINUOUS = True` (or pass `--watch`) to keep polling every `POLL_SECONDS` and upload projects as
  they reach `done` — a worker to run alongside a still-processing dataset. `--once` forces a single
  pass. Ctrl-C stops the watch loop.
- Because it selects off `done`, you never upload a half-processed project.
- Its resume file only lists projects whose load finished. A project with a `.project_id` but missing
  from the resume file was interrupted mid-upload.

## `zip_projects.py` — archive and reclaim disk

Zips every `done` project (optionally minus an exclude list) into `ARCHIVE_DIR`, then deletes the
original directory.

- ⚠️ **Destructive**: it `rmtree`s each project after zipping. Run it only once the data is safely
  uploaded/backed up. Set `DELETE_AFTER = False` for a dry run that creates the zips but keeps the
  originals — a good idea the first time.
- The zips keep each project's `.project_id`. Unzip the project back in place before re-uploading it,
  or the load creates a duplicate.
- `EXCLUDE_FILE` lets you hold back specific `rel_path`s (e.g. projects you're not done with).
- The archive directory is created automatically.

## How these relate to the core skill

`mwf ds load` and these scripts use the same `Dataset` object and state model the rest of the skill
relies on — they just read `get_dataframe(query_state=['done'], ...)` to find their work. That's the
payoff of the state model: once processing marks projects `done`, downstream steps consume that signal
without any separate bookkeeping. Keep the scripts' `DB_PATH` pointing at the same dataset file the
workflow updates.
