# IRB cluster (SLURM) specifics

Read this before writing or submitting any `sbatch` script on the IRB Barcelona cluster. The values
below were checked on the cluster (October 2026); re-check with the commands given when something
looks different. Paths use placeholders: `<conda>` is the user's conda installation (e.g.
`~/miniforge3`), `<env>` the environment with the workflow installed (e.g. `mwf_env`), and `<group>`
the user's SLURM group/account. Ask the user, or check their existing job scripts, instead of guessing.

## Access and environment

- Log in through the user's ssh alias for the cluster (often `ssh irb`). Large outputs go to
  `/scratch/<group>/$USER` (often linked as `~/scratch`), a BeeGFS filesystem shared by everyone and
  often >90% full: check `df -h /scratch` before a large run.
- **Scratch deletes old files.** A cleanup policy removes files on scratch by their timestamp. Copy
  tools that keep the source modification time make fresh downloads look old, so they get deleted
  during the run. This applies to `rclone copy`, `cp -p`, `rsync -a` and `tar`/`unzip`. **Update
  the timestamps after every copy into scratch**, in the same job as the copy:
  ```bash
  find <dataset>/<projects> -type f -exec touch {} +
  ```
  Symptom: projects that were downloaded in full later fail with `InputError: Input topology ... is
  nowhere to be found`, and their trajectory dirs are empty. In MDBind this hit 569 of 2157 projects
  in one array run. Once deleted, the files must be downloaded again.
  It also hits datasets that are already processed: in NMRlipids (FAIRMD BilayerData, 882 projects)
  every raw input, every `README.yaml` and the files under the cloned repo's `.git` were deleted a
  few weeks after processing, while the newer workflow outputs survived. The directories stay, so the
  dataset looks intact until a re-run (`mwf run -ow`) fails. Before planning a re-processing, check
  that the files named in each `inputs.yaml` still exist. Anything that must outlive the scratch
  cleanup (the dataset file, the `.project_id` upload traces, the dataset README, the helper
  scripts) belongs outside scratch or needs a backup there.
- `module` is **not available** on the compute nodes (non-login shell), so `module load anaconda3`
  fails. Activate conda in every job script with:
  ```bash
  source <conda>/etc/profile.d/conda.sh
  conda activate <env>
  ```
  For non-interactive commands on the login node, either source the same file or use `bash -lc "..."`.
- If the environment has an editable install of the workflow (`pip show mddb_workflow`, or
  `python -c "import mddb_workflow; print(mddb_workflow.__file__)"`), jobs run whatever is checked out
  in that repository. Changing it affects every running and future job; that is the user's call.
- Do not run heavy work on the login node. For a one-project test use `srun` (it blocks and streams the
  output), and for anything long (pre-flight, `add`, `inputs` over thousands of projects) an `sbatch` job:
  ```bash
  cd <project> && srun -p <lab partition> -c 2 --mem=8G -t 01:00:00 bash -c \
    "source <conda>/etc/profile.d/conda.sh && conda activate <env> && mwf run ..."
  ```

## Partitions

| Partition | Nodes | PriorityTier | Preemption | Allowed QoS |
|-----------|-------|--------------|------------|-------------|
| `mmb_cpu_zen3` | irbccn[10-12,14-15], 5 nodes, 280 CPUs (zen3) | 2 | OFF | `short`, `medium`, `long`, `interactive` |
| `mmb_cpu_sphr` | irbccn[17-22,50], 7 nodes, 448 CPUs (Sapphire Rapids) | 2 | OFF | `short`, `medium`, `long`, `interactive` |
| `spot_cpu` | irbccn[10-12,14-53], 43 nodes, 5392 CPUs (zen3, sph, zen4) | 1 | REQUEUE | `spot` only |

- **Lab partitions first.** Partitions are owned by labs (`AllowGroups`/`AllowAccounts` in
  `scontrol show partition`), and members get priority and no preemption there. The table shows the MMB
  ones (`mmb_cpu_zen3`, `mmb_cpu_sphr`: groups `orozco,its`, accounts `mmb,its`); for another lab, list
  the partitions it can use with `sacctmgr -n -P show assoc user=$USER format=Account,Partition,QOS` and
  `sinfo -s`. Use the lab partitions by default rather than the cluster default (`irb_cpu_zen4`).
  Partitions with the same QoS can be combined in one job, so tasks start on whichever has free CPUs:
  ```bash
  #SBATCH --partition=mmb_cpu_zen3,mmb_cpu_sphr      # MMB members
  ```
- **QoS on lab partitions is assigned automatically.** Array tasks with `--time=02:00:00` ran with QoS `short`
  (MaxWall 04:00:00), and the setup job on spot with `spot`. Keeping per-task time under 4 h stays
  within `short`; check what a job got with `sacct -j <id> -X -o JobID,Partition,QOS`.
- **`spot_cpu` runs on idle nodes of any lab** at lower priority. When the owners need their nodes,
  the job is preempted and requeued (cluster `PreemptMode=REQUEUE`), and the batch script **restarts
  from the beginning**. Make spot jobs safe to restart:
  ```bash
  #SBATCH --partition=spot_cpu
  #SBATCH --requeue
  #SBATCH --open-mode=append      # keep the log of the run that was preempted
  ```
  The `done`-state guard in `assets/job_array.sh` plus the workflow's own task cache make a restarted
  `mwf run` cheap. Setup scripts must also be idempotent: `mwf dataset add` skips UUIDs already in the
  dataset, `mwf dataset inputs` without `-o` skips existing files, and any summary file a script writes
  should be overwritten rather than appended to. A specific architecture can be requested with
  `--prefer=zen4` or `--nodelist=...`.
- **Spot cannot be a fallback inside a lab-partition job.** `spot_cpu` only allows QoS `spot` and the
  lab partitions do not allow it, so any `--partition` list mixing them is rejected with
  `allocation failure: Invalid qos specification`, whatever `--qos` is forced. To use both, submit
  separate jobs (for example a lab-partition array and a spot array working on separate lists).
- **Check a partition/QoS combination before submitting** with `sbatch --test-only ...`: it reports
  where and when the job would start, or the error.
- **The partition cannot be changed after submission** (`scontrol update ... Partition=` is rejected).
  Cancel pending jobs and resubmit.
- Inspect the current setup with:
  ```bash
  scontrol show partition <name> | grep -oE "PriorityTier=\S+|PreemptMode=\S+|AllowQos=\S+|Nodes=\S+"
  sacctmgr -n -P show qos short,spot format=Name,Priority,MaxWall
  sinfo -p mmb_cpu_zen3,mmb_cpu_sphr,spot_cpu -o "%P %a %D %C"     # CPUs allocated/idle/other/total
  sinfo --format="%13P %25N %5D %4c %b" -p spot_cpu                # node features (zen3/sph/zen4)
  ```

## Array limits

- `MaxArraySize = 6001`, so one array holds at most 6000 tasks (indices 0–6000). For larger datasets
  split the project list into parts and chain one array per part with `--dependency=afterany:<prev>`,
  so a single concurrency cap (`%K`) applies to the whole run:
  ```bash
  split -l 6000 -d --additional-suffix=.txt lists/all.txt lists/all_part
  DEP=""
  for part in lists/all_part*.txt; do
      N=$(wc -l < "$part")
      JOB=$(sbatch --parsable $DEP --array=0-$((N - 1))%40 --export=ALL,LIST=$part job_array.sh)
      DEP="--dependency=afterany:$JOB"
  done
  ```
  `MaxJobCount = 50000`. Check both with `scontrol show config | grep -E "MaxArraySize|MaxJobCount"`.
- Index arrays into a **list file** (one project dir per line, `sed -n "$((SLURM_ARRAY_TASK_ID + 1))p"`)
  instead of a live `ls` glob, so the mapping cannot shift between submissions.

## Submitting

- **Submit from the dataset root.** `--output=array_logs/...` is relative to the submit directory; if
  that directory does not exist the job fails immediately and writes no log at all (a job submitted
  from `utils/` died in 3 s this way). `mkdir -p array_logs` first.
- Cancelling an array kills its running tasks; those projects show up as "Workflow process is no longer
  running" and belong in the next retry list.
- All compute nodes share one outbound IP. External services queried by the workflow (PubChem, the RCSB
  ligand page, NCBI BLAST) rate-limit per IP (HTTP 429), so the limit is on the **whole run's**
  concurrency, not per node. See "Workflow flags for large datasets" in `launching.md`.

## Local resources

- Local BLAST (`mwf run -lb`) uses the Swiss-Prot database at
  `<workflow repo>/mddb_workflow/resources/blastdb/swissprot.*` (BLAST+ must be in `<env>`). It is
  downloaded on first use; make sure it exists before launching an array so tasks do not race to
  download it (see `launching.md`).
