#!/bin/bash
#SBATCH --partition=mmb_cpu_zen3,mmb_cpu_sphr
#SBATCH --job-name=ds_download
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00

# Download selected files from an S3/MinIO bucket with rclone into the dataset's project folder.
# Adaptable example (from the MDBind JSC batch: 23,771 files, 707 GiB, 35 min at ~300 MiB/s).
# See references/data-lifecycle.md, "S3/MinIO buckets with rclone".
#
# Needs: the rclone remote configured on the machine running this (~/.config/rclone/rclone.conf),
# and FILES_LIST = one path per line, relative to REMOTE (e.g. <id>/<id>.parm7), built from a
# saved listing so the bucket is not listed again.
# Launch from the dataset root: sbatch utils/download_rclone.sh
# Re-running is safe: files already present with the same size are skipped.

# --- CONFIG (adapt per dataset) ---------------------------------------------------------------
BASE="/home/<user>/scratch/<dataset>"       # dataset root (absolute: SLURM copies the script)
REMOTE="<remote>:<bucket>/<path>"           # rclone remote + path to the projects
DEST="$BASE/projects"                       # where <id>/... lands (keep the bucket layout)
FILES_LIST="$BASE/utils/download_files.txt" # the needed files only
EXCLUDED_DIR="$BASE/excluded"               # projects moved out of DEST must not come back
# ----------------------------------------------------------------------------------------------

TAG=${SLURM_JOB_ID:-local}
FILES=$BASE/utils/download_files_${TAG}.txt

# Drop projects that were excluded on purpose, otherwise they reappear in DEST and shift array indices
EXCL=$( [ -d "$EXCLUDED_DIR" ] && ls "$EXCLUDED_DIR" | paste -sd'|' )
grep -vE "^(${EXCL:-__none__})/" "$FILES_LIST" > "$FILES"
echo "files to check: $(wc -l < "$FILES") (excluded: ${EXCL:-none})"

list_missing() {
    while read -r f; do [ -e "$DEST/$f" ] || echo "$f"; done < "$FILES" > "$1"
    echo "$(date '+%F %T') missing files: $(wc -l < "$1") -> $1"
}
list_missing "$BASE/utils/download_missing_${TAG}_before.txt"

# --no-traverse: don't list the bucket, fetch exactly the files in the list.
# --local-no-set-modtime + --size-only: give files the current time (scratch purges by timestamp,
# and rclone would otherwise keep the old remote mtime) and stop comparing mtimes on re-runs.
rclone copy "$REMOTE" "$DEST" \
    --files-from "$FILES" \
    --no-traverse \
    --size-only --local-no-set-modtime \
    --transfers 16 --checkers 32 \
    --stats 5m --stats-one-line \
    --log-level INFO --log-file "$BASE/utils/download.log"

# Refresh timestamps of everything already on disk (older downloads and workflow outputs too)
find "$DEST" -type f -exec touch {} +
echo "$(date '+%F %T') timestamps refreshed"
list_missing "$BASE/utils/download_missing_${TAG}_after.txt"
