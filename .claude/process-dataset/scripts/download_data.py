"""Download a dataset's raw project archives and unpack them into project directories.

Adaptable example (based on the ATLAS setup). It builds a link list, downloads with aria2c
(resumable, segmented), then unzips each archive into its project folder. The SOURCE-SPECIFIC parts
are the CONFIG block below and `unzip_project` (which knows this source's archive layout) — adapt those
for a different provider; the aria2c mechanics are reusable as-is.

Requires the `aria2c` binary. Usage:  python download_data.py
"""
from pathlib import Path
import subprocess
import zipfile
import shutil

# --- CONFIG (adapt per data source) ---------------------------------------------------------------
API_BASE = "https://www.dsimb.inserm.fr/ATLAS/api"       # base URL of the download API
ID_LIST_PATH = Path(__file__).parent / "id_list.txt"      # one entry ID per line
LINKS_FILE = Path(__file__).parent / "links.txt"          # generated aria2c input file
DEST_DIR = Path(__file__).parent.parent / "projects"      # where project dirs are created
# --------------------------------------------------------------------------------------------------


def read_id_list(txt_file: Path) -> list[str]:
    """Read non-empty entry IDs from a text file."""
    with open(txt_file) as f:
        return [line.strip() for line in f if line.strip()]


def download_with_aria2c(retries: int = 5, segments: int = 5, verify_ssl: bool = True) -> None:
    """Download every URL in LINKS_FILE with aria2c, resuming partial files."""
    cmd = [
        "aria2c",
        "--continue=true",
        "-i", LINKS_FILE.as_posix(),
        f"--max-tries={retries}",
        f"--max-connection-per-server={segments}",
        "--remote-time=true",
        "--content-disposition=false",
        "-d", DEST_DIR.name,
    ]
    if not verify_ssl:
        cmd.append("--check-certificate=false")
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"aria2c failed with code {result.returncode}")


def unzip_project(zip_path: Path) -> None:
    """Unpack one downloaded archive into its project directory (source-specific layout)."""
    print(f"Unzipping {zip_path}...")
    proj_dir = DEST_DIR / zip_path.stem[:6]
    proj_dir.mkdir(parents=True, exist_ok=True)
    extract_path = proj_dir / "_extract"
    extract_path.mkdir(parents=True, exist_ok=True)

    try:
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(extract_path)
    except zipfile.BadZipFile:
        shutil.rmtree(extract_path)
        zip_path.unlink(missing_ok=True)
        print("Error extracting. Removing corrupted zip.")
        return

    # This source nests useful files under a "protein/" subfolder; fall back to the archive root.
    protein_dir = extract_path / "protein"
    source_dir = protein_dir if protein_dir.exists() else extract_path

    for item in source_dir.iterdir():
        target = proj_dir / item.name
        if target.exists():
            if target.is_dir() and item.is_dir():
                for nested_item in item.iterdir():
                    shutil.move(str(nested_item), str(target / nested_item.name))
                item.rmdir()
            else:
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()
                shutil.move(str(item), str(target))
        else:
            shutil.move(str(item), str(target))

    shutil.rmtree(extract_path)
    zip_path.unlink(missing_ok=True)


def main():
    ids = read_id_list(ID_LIST_PATH)
    urls = []
    for entry in ids:
        dest_path = DEST_DIR / entry
        # Skip entries already downloaded (folder present, or archive already fetched)
        if dest_path.exists() or dest_path.with_stem(entry + '_protein.zip').exists():
            continue
        urls.append(f"{API_BASE}/ATLAS/protein/{entry}")

    LINKS_FILE.write_text("\n".join(urls))
    print(f"Prepared {len(urls)} / {len(ids)} download links in {LINKS_FILE}.")
    download_with_aria2c()

    zip_list = sorted(DEST_DIR.glob("*.zip"))
    print(f"Unzipping {len(zip_list)} projects.")
    for file in zip_list:
        unzip_project(file)


if __name__ == '__main__':
    main()

# Tip — find project folders that ended up empty (only a cache file):
#   find projects -type d -exec bash -c '[[ $(ls -A "$1" | wc -l) -eq 1 && -f "$1/.mwf_cache.json" ]]' _ {} \; -print
