"""Archive finished projects into per-project zips and reclaim disk space, driven by dataset state.

For every project in state 'done' (optionally minus an exclude list), create a zip under ARCHIVE_DIR
and then delete the original project directory. This is the "archive" step once a project has been
processed and uploaded — the dataset's `done` state is the queue.

WARNING: this DELETES each project directory after zipping it (irreversible). Run it only after the
data is safely uploaded/backed up, and consider a dry run first (set DELETE_AFTER = False).

Usage:  python zip_projects.py
"""
import subprocess
import shutil
from pathlib import Path

from mddb_workflow.core.dataset import Dataset

# --- CONFIG ---------------------------------------------------------------------------------------
DB_PATH = "../mwf_ds.db"                    # dataset file
SUBMIT_DIR = Path("..").resolve()           # dataset root (parent of the project dirs)
ARCHIVE_DIR = SUBMIT_DIR / "zipped_completed"
EXCLUDE_FILE = None                         # optional path to a file of rel_paths to skip (or None)
DELETE_AFTER = True                         # set False for a dry run that zips but keeps originals
# --------------------------------------------------------------------------------------------------

ds = Dataset(DB_PATH)
ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)

done_paths = ds.get_dataframe(query_scope='projects', query_state='done')['rel_path'].to_list()

exclude = set()
if EXCLUDE_FILE and Path(EXCLUDE_FILE).exists():
    exclude = {line.strip() for line in Path(EXCLUDE_FILE).read_text().splitlines() if line.strip()}

to_zip = sorted(p for p in done_paths if p not in exclude)
print(f"Total done: {len(done_paths)} | skipping (excluded): {len(done_paths) - len(to_zip)} | to zip: {len(to_zip)}")

for rel_path in to_zip:
    project_dir = SUBMIT_DIR / rel_path
    if not project_dir.exists():
        print(f"Warning: {project_dir} does not exist, skipping.")
        continue

    zip_path = ARCHIVE_DIR / f"{rel_path.replace('/', '_')}.zip"
    print(f"Zipping {rel_path} -> {zip_path.name}")
    result = subprocess.run(
        ["zip", "-r", "-q", str(zip_path), rel_path],
        cwd=SUBMIT_DIR, text=True, capture_output=True,
    )
    if result.returncode == 0:
        if DELETE_AFTER:
            shutil.rmtree(project_dir)
            print(f"Archived and removed {rel_path}")
        else:
            print(f"Archived {rel_path} (kept original; DELETE_AFTER is False)")
    else:
        print(f"Error zipping {rel_path}: {result.stderr}")

print(f"\nDone. Archives in {ARCHIVE_DIR}")
