# Data lifecycle: download, upload, clean up

Processing is the middle of a longer pipeline. Uploading is built into the CLI (`mwf ds load`);
downloading, measuring and archiving use helper scripts (standalone executables in `scripts/`). Upload
and archive are driven by the dataset's `done` state, which makes the dataset the natural work queue: a
project becomes eligible to publish/archive exactly when the workflow marks it `done`.

```
  download  ─►  add + inputs + run   ─►  mwf ds load    ─►  clean up
 (acquire)      (this skill's core)      (publish done)     (measure, classify, archive on RES, delete)
```

The scripts are **adaptable examples**, not turnkey tools — most have a CONFIG block at the top with
the dataset-specific paths to set. Upload and cleanup are outward-facing or irreversible; read the
warnings.

Contents:
- `download_data.py` — acquire raw projects
- S3/MinIO buckets with rclone
- `mwf ds load` — publish finished projects; the upload record (`.loader_register.json`), uploading
  from another host
- `upload_projects.py` — legacy upload script
- Clean up after the upload — checks, inventory, classes, public or shared source, deleting;
  archiving on RES (MareNostrum 5); making the zips with `zip_projects.py`
- How these relate to the core skill

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
mwf ds ls mwf_ds.sql -st done -n 0                              # what would be uploaded
dataset_inventory.py . --uploads                                # which ones are uploaded already
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

### The upload record: `.loader_register.json`

The loader leaves `.loader_register.json` in each project directory it loads: a list with one entry
per load (date, loader version, command line, project id, accession). Older loaders left `.project_id`,
with the id alone; the current one still reads it when there is no register, and deletes it once it
writes the register. When the loader starts on a project, it decides which database entry to use, in
this order:

1. an accession forced on the command line or in the inputs file (`accession:`);
2. the project id of the last register entry (or of `.project_id`), if that project still exists in
   the database;
3. otherwise, a new project with a new accession.

So re-uploading a re-processed project updates the same entry, **as long as the record is still
there**. Consequences:

- **Never delete `.loader_register.json` or `.project_id` files**, and don't delete or recreate
  project directories that hold one. If the file is lost, the next load creates a duplicate project.
  On a scratch filesystem with an age-based cleanup, back them up (with their relative paths) outside
  scratch, or write the accession into each `inputs.yaml` (`mwf ds inputs -s accession <A0XXX> -p
  <project>`) so it no longer depends on the file.
- If the id points to a project that no longer exists in the database, the loader warns, removes the
  record and creates a new project.
- The record is written when the load **starts**, so a project with one may still be a partial
  upload. Compare with the database (the MDDB tools or the web client) before assuming it is
  complete, and re-run `load` on it.
- **The records stay on the copy that ran the load.** When a dataset is uploaded from another host
  (a copy made for the upload, the cluster), copy the records back to the main copy afterwards, or the
  next load from there creates every project again. Copy only those files, with a dry run first
  (`-n --itemize-changes`; IDRome: 27,998 new, 60 updated, 6 s):
  ```bash
  rsync -rtp -O --include='*/' --include='.loader_register.json' --exclude='*' --prune-empty-dirs \
      <host>:<root>/ <root>/
  ```
  Records that exist on both sides should only have gained entries with the same project id; check
  that before overwriting.
- `mwf ds ls` only looks for `.project_id` (`Loaded: yes`), so it says `no` for projects uploaded by
  the current loader. `scripts/dataset_inventory.py <root> --uploads` checks the records instead (see
  "Clean up after the upload").

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
- Its resume file only lists projects whose load finished. A project with an upload record but
  missing from the resume file was interrupted mid-upload.

## Clean up after the upload

After the upload, most of a dataset's disk is extra copies, files the source still has, and outputs
the node already holds. Deleting them can't be undone, so the cleanup is a step of its own, with its
own checks. Worked example: IDRome (28,058 projects) took 134 GB on the workstation, plus a 73 GB
copy made on another host for the upload; 61 GB of it was the authors' archive and files the workflow
never read.

**1. Check first.** Delete nothing until all of these hold:
- **The upload is complete.** `scripts/dataset_inventory.py <root> --uploads` reports no problem:
  every project (a folder with an `inputs.yaml`) has an upload record, no accession is shared by two
  projects, and every record has at least one complete load. A load whose recorded `command` has
  `--include`, `--exclude`, `--skip-trajectories`, `--skip-files`, `--skip-analyses` or
  `--md-directories` uploaded only part of the project. No load failed either (the load logs, or the
  `mwf ds load` exit codes), and a few projects read back from the node's API look right. The node's
  project count (`<node>/api/rest/v1/projects?limit=1`, field `filteredCount`) is only a rough check:
  it may count only published projects, and other uploads change it.
- The upload records are on the copy future loads will run from (see "The upload record" above).
- The projects are published, if that was the plan, and the data owners have no corrections pending.
- The user has said which hosts hold copies, which one is the main copy (the one future loads run
  from), and where the data came from (see "Public or shared source" below).

**2. Measure every copy.** `scripts/dataset_inventory.py <root>` walks the dataset without following
links and prints the size per file name, largest first. Digit runs count as one name
(`mwf_3.out` and `mwf_12.out` are `mwf_N.out`), and `--by ext` groups by extension instead. It read
IDRome's 1.6 million files in 13 s on a local disk. On another host, run it without copying it:
`ssh <host> python3 - <root> < dataset_inventory.py`.

**3. Sort the files into classes.** Write the classes in a file next to the dataset
(`utils/cleanup_classes.yaml`, which also records the decision; YAML needs PyYAML, which `mwf_env`
has, or write JSON) and re-run with `--classes` until nothing large is left in `unclassified`:
- A pattern without `/` matches the file name at any depth.
- A pattern with `/` matches the path relative to the root, and its `*` also crosses `/`:
  `IDRome_v4/*/mwf_*.out` is every `mwf_*.out` below `IDRome_v4/`. Use such paths to keep a class
  inside the project folders.
- A pattern ending in `/` is a whole folder (`array_logs/`).
- Patterns match files, not dataset states: a class takes the files of every project, `done` or not.
  If some projects are not `done`, leave them out by path.

The script never lists its protected files: the upload records (`.loader_register.json`,
`.project_id`), `.mwf_cache.json`, `.register.json`, `inputs.yaml`, `metadata.json`, `topology.json`,
any `*.db` or `*.sql`, `README.md`, and everything in `utils/`, `run_logs/` and `array_logs/` at the
root. Add more with a `protected:` list in the classes file (other log folders, notes). A class that
matches a protected file is an error, and `--list` prints nothing until it is fixed. The script also
reports files that two classes match and patterns that match nothing. Watch the symlinks: the workflow
may link `replica_1/trajectory.xtc` to the raw trajectory, and `--list` warns about links that would
break. Put such links in the same class as their target, or keep the target.

**Public or shared source.** Where the data came from decides what happens to the files we received
(the raw inputs, and the source files the workflow never read):
- **Public and lasting:** a public URL, a repository with a DOI (Zenodo, Figshare), a public
  database. Delete them after the upload: they can be downloaded again, so the README only has to say
  where and how (URL or DOI, download command). IDRome (Zenodo) is this case; nothing of it is
  archived.
- **Shared directly with us:** a hard drive, an S3/MinIO bucket (MDBind), a folder shared with us on
  a server, a private or expiring link. Archive them on RES ("Archiving on RES" below) before
  deleting the local copy: we may be the only ones who still have them.

When it is not clear which case it is, ask the user; never assume a share will stay online.

| Class | Typical files | Default |
|---|---|---|
| bookkeeping | the protected files above, and anything else the dataset's README relies on | always keep |
| extra copies | a copy made on another host for the upload, a scratch copy on the cluster | delete, after copying back their logs and upload records |
| unused source files | the source archive once extracted, files the workflow never read (IDRome: `IDRome.zip`, `emap.csv`, `*.npy`) | public source: delete. Shared directly: archive on RES, then delete |
| raw inputs | the topologies and trajectories named in `inputs.yaml` | public source: delete; the workflow downloads nothing, so download them again (the README has the command) before any re-run. Shared directly: archive on RES, then delete |
| caches and logs | `.*_offsets.npz`, and the projects' `mwf_*.out` and `mwf_*.err` | delete |
| workflow outputs | `mda.*.json` analyses, `mdt.*.xtc`, structures, screenshots | delete once every load was complete and no more corrections are expected: they are on the node, and a re-upload should re-run the whole workflow anyway (RES README) |

```yaml
# utils/cleanup_classes.yaml (IDRome)
caches: ['.*_offsets.npz', 'IDRome_v4/*/mwf_*.out', 'IDRome_v4/*/mwf_*.err']
source_unused: ['IDRome.zip', 'emap.csv', '*.npy', 'restart.gsd', 'atommass.dat', 'analysis.csv', 'top_CA.pdb']
raw_inputs: ['IDRome_v4/*/traj.xtc', 'IDRome_v4/*/top.pdb', 'IDRome_v4/*/replica_1/trajectory.xtc']
outputs: ['mda.*.json', 'mdt.*.xtc', 'structure.cif', 'average.pdb', 'first_frame.pdb', 'mdf.screenshot.jpg']
```

**4. Decide each class with the user.** Show the class table with the size per host, whether the
source is public or was shared directly, and the default for each class. Write the decision into the
README.

**5. Archive on RES what must be kept.** For the classes from a source shared directly with us, make
the zips, copy them to RES and check them there ("Archiving on RES" below). Only then is their local
copy deleted, in step 6.

**6. Delete one class at a time.**
1. First, while the main copy still has everything, copy back from each extra copy what only it has
   (its logs and upload records; see "A whole extra copy" below), so the backup holds them too.
2. Once, before the first class, archive the protected files (the bookkeeping) on RES. Zip them
   outside the dataset folder, with the log outside too: a log written into the dataset while it is
   zipped makes zip warn that a file changed.
   ```bash
   cd <root> && set -o pipefail && dataset_inventory.py . --classes utils/cleanup_classes.yaml \
       --list protected --relative | zip -T -@ <outside>/<dataset>_bookkeeping_<date>.zip > <outside>/zip.log 2>&1 \
     && grep -c 'zip warning' <outside>/zip.log; zipinfo -1 <outside>/<dataset>_bookkeeping_<date>.zip | wc -l
   ```
   The zip must pass its test with no warnings, and the count must match the `protected` row of the
   class table. Then copy it to RES and check it as in "Archiving on RES" below (the README's rsync
   command, MD5 on both sides), into the folder the archive's README gives for uploaded data, with a
   short note next to it: what the zip holds, where the data itself is (DOI, node), its MD5. The
   bookkeeping compresses well (IDRome: 1.5 GB of JSON in 366 MB), so the zip is usually under RES's
   1 GB rule: ask the user whether to store it uncompressed (`zip -0`) instead. Delete the local zip
   once the RES copy is checked. The README's cleanup entry is written after the backup, so the note
   next to the zip is what explains it.
3. List the class and check what the script prints (on stderr) against the class table: path count,
   size, patterns that match nothing, and any broken-link warning. Look at a few paths.
   ```bash
   dataset_inventory.py <root> --classes utils/cleanup_classes.yaml --list caches -0 > caches.list
   ```
4. Confirm with the user: the class, its count and size, and that it can't be undone. Then delete
   exactly that list: `xargs -0 rm -- < caches.list`. Never `rm -r` project directories: they hold
   the upload records.
5. Re-run the inventory: the class is empty, and `protected` has the same count as before.

**A whole extra copy** goes with one `rm -rf` of its root instead:
- First copy back its logs and upload records (item 1 above), and check that it holds nothing the main
  copy lacks: `rsync -rtp -O -n --itemize-changes <host>:<root>/ <root>/` should list only files you
  chose to leave behind. Without `-l`, rsync reports every symlink as a skipped "non-regular file". If
  the classes are already deleted here, the dry run lists them too: sort its paths with the classes
  file and look only at the rest (IDRome: 193 of 1.15 million). Compare the dataset files by their
  rows (Python `sqlite3` `iterdump()` on read-only copies); in IDRome the upload had changed only the
  remembered `loader_directory` setting.
- Check that it is not the main copy: future loads run from the copy that keeps the records.
- Confirm with the user the host, the absolute path and its size (`du -sh`), and that it can't be
  undone. Then delete it with the absolute path written out in full (`ssh <host> rm -rf <path>`),
  never through a variable that could be empty.

**7. Record it in the README:** the date, what was deleted on which host and its size, and where each
class can be fetched again (the public URL or DOI, the node, the RES archive path).

### Archiving on RES (MareNostrum 5)

Data that must be kept (a source shared directly with us, see "Public or shared source") and the
bookkeeping backup of step 6 go to the group's RES archive on MareNostrum 5. Its path and the transfer node that reaches it are not written
in this skill: ask the user. Before writing anything there, read the `README` at the root of the
archive. It says which folder each kind of data goes to, how to record what was archived, and the
rules below, which may have changed since they were copied here (2026-11):

- **Only files of 1 GB or more.** The quota is 10% SSD and 90% tape. A file over 1 GB that is not
  read or written for a month moves to tape; small files never do, and they fill the SSD quota (it
  happened). Group projects into zip or tgz files of at least 1 GB (one large project, one system or
  a range of IDs per file), never loose project folders.
- **Check every zip.** Make it with `zip -rT` (tests the archive), without `-q`, which also hides the
  warnings, and read its output: `zip warning: file size changed while zipping` means a file is
  incomplete in the zip, so make it again. `zip_projects.py` does both.
- **Don't keep the source dates when copying** (no `rsync -a`), or the new files look a month old and
  go to tape at once. The README's command:
  ```bash
  rsync -vrlpgoDP --size-only --append <archives> <transfer node>:<archive folder>/
  ```
- **Verify before deleting the local copy:** the same size and MD5 on both sides (`md5sum` on the
  transfer node right after the copy, while the files are still on SSD).
- **Record it** as the archive's README asks (which accessions each archive holds, with notes on the
  data), and the archive names and MD5 in the dataset's README.

### Making the zips with `zip_projects.py`

`zip_projects.py` zips every `done` project (optionally minus an exclude list) into
`ARCHIVE_DIR/<rel_path>.zip`. Each zip is made from scratch and tested with `zip -T`; when zip warns
that a file changed while it was read, the zip is deleted, the project is kept, and the script exits
with 1 and names it. It makes one zip per project, which suits RES only when every project is over
1 GB; otherwise group smaller projects into larger zips (by ID prefix, by batch).

- By default it keeps the projects (`DELETE_AFTER = False`): in the cleanup their files are deleted
  later (step 6), through a class list, once the zips are checked on RES.
- ⚠️ `DELETE_AFTER = True` removes each whole project folder right after its zip passes, upload
  records and bookkeeping included. Use it only when the zips themselves are the final copy.
- The zips keep each project's upload record. Unzip the project back in place before re-uploading it,
  or the load creates a duplicate.
- `EXCLUDE_FILE` lets you hold back specific `rel_path`s (e.g. projects you're not done with).
- The archive directory is created automatically.

## How these relate to the core skill

`mwf ds load`, `upload_projects.py` and `zip_projects.py` use the same `Dataset` object and state
model the rest of the skill relies on — they read `get_dataframe(query_state=['done'], ...)` to find
their work. That's the payoff of the state model: once processing marks projects `done`, downstream
steps consume that signal without any separate bookkeeping. Keep the scripts' `DB_PATH` pointing at
the same dataset file the workflow updates. `dataset_inventory.py` is the exception: it reads files,
not the dataset, so it doesn't know which projects are `done`.
