"""Zip every finished project of a dataset, one zip per project, driven by dataset state.

For every project in state 'done' (optionally minus an exclude list), make ARCHIVE_DIR/<rel_path>.zip
(the archive folder mirrors the project paths, so two projects never get the same zip name). Each zip
is made from scratch, tested (`zip -T`), and rejected when zip warns that a file changed while it was
read (that file would be incomplete in the zip). A rejected zip is deleted, so a re-run makes it again.

By default the projects are kept (DELETE_AFTER = False). In step 7 of the skill (clean up after the
upload) this script only makes the zips: they are copied to the archive and checked there (MD5) first,
and the project files are deleted afterwards through a reviewed class list (scripts/dataset_inventory.py),
never by this script. DELETE_AFTER = True removes each whole project folder right after its zip passes,
upload records and bookkeeping included; use it only when the zips themselves are the final copy.

For the RES archive on MareNostrum 5 the zips must be 1 GB or more: group small projects into larger
zips instead of one zip per project (see "Archiving on RES" in references/data-lifecycle.md).

Usage:  python zip_projects.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

from mddb_workflow.core.dataset import Dataset

# --- CONFIG ---------------------------------------------------------------------------------------
DB_PATH = "../mwf_ds.db"                    # dataset file
SUBMIT_DIR = Path("..").resolve()           # dataset root (parent of the project dirs)
ARCHIVE_DIR = SUBMIT_DIR / "zipped_completed"
EXCLUDE_FILE = None                         # optional path to a file of rel_paths to skip (or None)
DELETE_AFTER = False                        # True deletes each project folder once its zip passes
# --------------------------------------------------------------------------------------------------

ds = Dataset(DB_PATH)
ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)

done_paths = ds.get_dataframe(query_scope='projects', query_state='done')['rel_path'].to_list()

exclude = set()
if EXCLUDE_FILE and Path(EXCLUDE_FILE).exists():
    exclude = {line.strip() for line in Path(EXCLUDE_FILE).read_text().splitlines() if line.strip()}

to_zip = sorted(p for p in done_paths if p not in exclude)
print(f"Total done: {len(done_paths)} | skipping (excluded): {len(done_paths) - len(to_zip)} | to zip: {len(to_zip)}")

failed = []
for rel_path in to_zip:
    project_dir = SUBMIT_DIR / rel_path
    if not project_dir.exists():
        print(f"Warning: {project_dir} does not exist, skipping.")
        continue

    zip_path = ARCHIVE_DIR / f"{rel_path}.zip"
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    # Always start from an empty zip: zip -r would otherwise update a zip left by an earlier run
    zip_path.unlink(missing_ok=True)
    print(f"Zipping {rel_path} -> {zip_path.relative_to(ARCHIVE_DIR)}")
    # Not -q: quiet mode also hides the "file size changed while zipping" warning checked below
    result = subprocess.run(
        ["zip", "-r", "-T", str(zip_path), rel_path],
        cwd=SUBMIT_DIR, text=True, capture_output=True,
    )
    output = result.stdout + result.stderr
    # zip only warns (exit 0) when a file changes while it is read, and that file is incomplete in the zip
    changed = 'file size changed while zipping' in output
    if result.returncode == 0 and not changed:
        if DELETE_AFTER:
            shutil.rmtree(project_dir)
            print(f"Archived and removed {rel_path}")
        else:
            print(f"Archived {rel_path} (kept the original)")
    else:
        zip_path.unlink(missing_ok=True)
        failed.append(rel_path)
        reason = 'a file changed while zipping' if changed else f'exit {result.returncode}'
        warnings = '\n'.join(line for line in output.splitlines() if 'warning' in line or 'error' in line)
        print(f"Error zipping {rel_path} ({reason}): deleted the zip, kept the original\n{warnings}")

print(f"\nDone. Archives in {ARCHIVE_DIR}. {len(failed)} failed{': ' + ', '.join(failed) if failed else ''}")
sys.exit(1 if failed else 0)
