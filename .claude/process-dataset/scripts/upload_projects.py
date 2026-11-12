"""Upload finished projects to the MDDB database with the loader, driven by dataset state.

Finds projects in state 'done', uploads each one that hasn't been uploaded yet (via the MDDB `loader`),
and records it in a resume file so re-runs skip it. This is the natural "publish" step after the
workflow finishes: the dataset's `done` state is the queue.

Two modes (set CONTINUOUS below, or pass --watch / --once on the command line):
  - one-shot (default): upload every project currently 'done', then exit.
  - continuous (--watch): keep polling every POLL_SECONDS and upload projects as they reach 'done' —
    useful to run alongside a still-processing dataset so results publish as they finish. Ctrl-C stops.

WARNING: uploading publishes data to an external database — it is an outward-facing action. Point
LOADER / DB_PATH at the right places and test on one project before running the full loop.

Usage:  python upload_projects.py [--watch | --once]
"""
import os
import sys
import time
import subprocess
from pathlib import Path

from mddb_workflow.core.dataset import Dataset

# --- CONFIG ---------------------------------------------------------------------------------------
DB_PATH = "../mwf_ds.db"                       # dataset file
UPLOADED_FILE = "uploaded_projects.txt"        # resume log (one rel_path per line)
CONTINUOUS = False                             # True = keep polling; False = one pass then exit
POLL_SECONDS = 30                              # wait between checks when the queue is empty (continuous)
# Location of the MDDB loader entry point (first existing path wins)
LOADER_CANDIDATES = [
    os.path.expanduser("~/repo/loader/index.js"),
    os.path.expanduser("~/repo/MDDB/loader/index.js"),
]
# --------------------------------------------------------------------------------------------------

# Command-line override of the CONTINUOUS config flag
if "--watch" in sys.argv:
    CONTINUOUS = True
elif "--once" in sys.argv:
    CONTINUOUS = False

loader = next((p for p in LOADER_CANDIDATES if os.path.exists(p)), None)
if loader is None:
    raise FileNotFoundError(f"Loader not found in any of: {LOADER_CANDIDATES}")

ds = Dataset(DB_PATH)


def load_uploaded() -> set[str]:
    if os.path.exists(UPLOADED_FILE):
        return set(Path(UPLOADED_FILE).read_text().splitlines())
    return set()


def upload_pending() -> int:
    """Upload every not-yet-uploaded 'done' project once. Returns how many were uploaded."""
    uploaded = load_uploaded()
    to_upload = ds.get_dataframe(
        query_state=['done'], query_scope='project', sort_by='rel_path', asc=True,
    )['rel_path'].to_list()
    to_upload = [p for p in to_upload if p not in uploaded]

    if to_upload:
        print(f"Projects to upload: {len(to_upload)}")
    total = len(to_upload)
    for i, proj in enumerate(to_upload):
        print(f"Uploading project: {proj} ({i + 1}/{total})...")
        subprocess.run(f'node {loader} load . -c', cwd=ds.root_path / proj, shell=True, check=True)
        with open(UPLOADED_FILE, 'a') as f:
            f.write(f"{proj}\n")
        print(f"Uploaded project: {proj}. {total - i - 1} remaining...")
    return len(to_upload)


def main():
    if not CONTINUOUS:
        n = upload_pending()
        print(f"Done. Uploaded {n} project(s)." if n else "Nothing new to upload.")
        return

    # Continuous mode: keep polling for newly-'done' projects until interrupted.
    print(f"Continuous mode: polling every {POLL_SECONDS}s (Ctrl-C to stop).")
    while True:
        if upload_pending() == 0:
            print(f"No new 'done' projects. Checking again in {POLL_SECONDS}s...")
        time.sleep(POLL_SECONDS)


if __name__ == '__main__':
    main()
