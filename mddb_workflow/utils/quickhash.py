"""Sampled xxHash quick-check of the files, the same as in the biobb provenance RO-Crates.

The quickhash identifies a file without reading it whole: the XXH3-128 (seed 0) of the file size
(uint64 big-endian) followed by three samples of 32 MiB at 0, floor((size - n) / 2) and size - n,
or by the whole content of the files smaller than 3 samples (96 MiB).
It is an identity check, not an integrity check: changes outside the samples go undetected.

32 MiB samples hold complete frames of a 1 million atom trajectory in the usual formats (XTC ~3.8 MB
per frame, DCD, NetCDF or TRR ~12 MB), so a change to some atoms of every frame is detected,
and the structure and topology files (smaller than 96 MiB) are hashed whole.

WARNING: This is the same scheme as biobb_common.provenance.quickhash, keep both implementations in sync.
WARNING: Both are tested against the same golden vectors (test/data/input/quickhash_golden.json here).
"""
from os.path import getsize

import xxhash

QUICKHASH_ALGO = 'xxh3_128-sample'
QUICKHASH_VERSION = 1
QUICKHASH_SAMPLE_SIZE = 1 << 25  # 32 MiB
# Read size of the samples and of the files hashed whole
HASH_CHUNK_SIZE = 1 << 20


def quickhash(path: str, sample_size: int = QUICKHASH_SAMPLE_SIZE) -> str:
    """Get the quickhash of a file: the lowercase hex XXH3-128 (seed 0) of its size as uint64 big-endian,
    followed by its samples of sample_size bytes at 0, floor((size - n) / 2) and size - n (in that order),
    or by its whole content if it is smaller than 3 samples.
    Raise an OSError if the file is shorter than its size when it is read (e.g. truncated during the read).
    """
    if sample_size < 1: raise ValueError(f'The sample size must be positive: {sample_size}')
    size = getsize(path)
    hasher = xxhash.xxh3_128()
    hasher.update(size.to_bytes(8, 'big'))
    with open(path, 'rb') as file:
        if size < 3 * sample_size:
            _update(hasher, file, size)
        else:
            for offset in (0, (size - sample_size) // 2, size - sample_size):
                file.seek(offset)
                _update(hasher, file, sample_size)
    return hasher.hexdigest()


def _update(hasher, file, size: int):
    """Hash exactly size bytes from the current position of the file, read in chunks,
    so a file truncated during the read raises instead of hashing short data.
    """
    remaining = size
    while remaining:
        chunk = file.read(min(HASH_CHUNK_SIZE, remaining))
        if not chunk:
            raise OSError(f"Unexpected end of file {getattr(file, 'name', '')}: wanted {size} bytes, got {size - remaining}")
        hasher.update(chunk)
        remaining -= len(chunk)
