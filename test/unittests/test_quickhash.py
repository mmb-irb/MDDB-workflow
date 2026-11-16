"""Tests of the sampled xxHash quick-check of the files.

biobb_common implements the same scheme (biobb_common/provenance/quickhash.py) and checks it against
the same golden vectors (test/data/input/quickhash_golden.json here): keep both in sync.
"""
import json
import os

import pytest
import xxhash

from mddb_workflow.utils import quickhash as qh

pytestmark = pytest.mark.CI

GOLDEN_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'input', 'quickhash_golden.json')
with open(GOLDEN_PATH) as golden_file:
    GOLDEN = json.load(golden_file)
# Small sample size, so the files of the boundary tests are small
N = 16


def content(size: int) -> bytes:
    """Content of the golden vector files: bytes(range(256)) repeated and cut to size."""
    return (bytes(range(256)) * (size // 256 + 1))[:size]


def write_content(path, size: int):
    """Write the content of a golden vector file in 1 MiB blocks (a whole number of patterns), without holding it in memory."""
    block = bytes(range(256)) * 4096
    with open(path, 'wb') as file:
        for start in range(0, size, len(block)):
            file.write(block[:size - start])
    return path


def write(path, data: bytes) -> str:
    path.write_bytes(data)
    return str(path)


def flipped(data: bytes, offset: int) -> bytes:
    changed = bytearray(data)
    changed[offset] ^= 0xFF
    return bytes(changed)


def digest(data: bytes, hashed: bytes) -> str:
    """XXH64 of the size prefix of data followed by the hashed bytes: the hash of the sampled files."""
    return xxhash.xxh64_hexdigest(len(data).to_bytes(8, 'big') + hashed)


@pytest.mark.parametrize('vector', GOLDEN['vectors'], ids=lambda v: f"size{v['size']}-n{v['sample_size']}")
def test_golden_vectors(tmp_path, vector):
    path = write_content(tmp_path / 'file', vector['size'])
    assert qh.quickhash(str(path), vector['sample_size']) == vector['quickhash']
    # The largest ones are 128 MiB: not kept in the pytest temporary directories
    path.unlink()


def test_golden_scheme():
    """Same scheme and default sample size as the biobb implementation."""
    assert (GOLDEN['algo'], GOLDEN['version'], GOLDEN['default_sample_size']) == (qh.QUICKHASH_ALGO, qh.QUICKHASH_VERSION, qh.QUICKHASH_SAMPLE_SIZE)


def test_small_files(tmp_path):
    """Files smaller than 3 samples: the plain XXH64 of the file, without size prefix (the hash of the MDDB cksum)."""
    for size in (0, 1, 2 * N, 3 * N - 1):
        data = content(size)
        assert qh.quickhash(write(tmp_path / 'small', data), N) == xxhash.xxh64_hexdigest(data)
    # With the default sample size (32 MiB)
    assert qh.quickhash(write(tmp_path / 'small', content(100))) == xxhash.xxh64_hexdigest(content(100))
    assert qh.quickhash(write(tmp_path / 'empty', b'')) == 'ef46db3751d8e999'


def test_boundary_sizes(tmp_path):
    # 3n - 1: whole file, without size prefix
    data = content(3 * N - 1)
    assert qh.quickhash(write(tmp_path / 'small', data), N) == xxhash.xxh64_hexdigest(data)
    # 3n: the samples are exactly adjacent, they cover the whole file (after the size prefix)
    data = content(3 * N)
    assert qh.quickhash(write(tmp_path / 'adjacent', data), N) == digest(data, data)
    # 3n + 1: middle sample at floor((2n + 1) / 2) = n, the byte at 2n is not sampled
    data = content(3 * N + 1)
    assert qh.quickhash(write(tmp_path / 'odd', data), N) == digest(data, data[:2 * N] + data[2 * N + 1:])


def test_sensitivity_and_blind_spot(tmp_path):
    # Samples of a 167 bytes file: [0, 16), [75, 91) and [151, 167)
    data = content(10 * N + 7)
    original = qh.quickhash(write(tmp_path / 'file', data), N)
    for offset in (5, 80, 160):
        assert qh.quickhash(write(tmp_path / 'file', flipped(data, offset)), N) != original
    # Known blind spot: a change outside the samples is not detected
    assert qh.quickhash(write(tmp_path / 'file', flipped(data, 40)), N) == original


def test_size_sensitivity(tmp_path):
    """Files with identical samples but different sizes."""
    assert qh.quickhash(write(tmp_path / 'a', bytes(10 * N)), N) != qh.quickhash(write(tmp_path / 'b', bytes(10 * N + 1)), N)


def test_determinism(tmp_path):
    path = write(tmp_path / 'file', content(10 * N + 7))
    assert qh.quickhash(path, N) == qh.quickhash(path, N)


@pytest.mark.parametrize('size, claimed', [(10, 20), (40, 100)], ids=['whole-file', 'samples'])
def test_truncated_during_read(tmp_path, monkeypatch, size, claimed):
    """A file shorter than its size when read raises instead of hashing short data."""
    path = write(tmp_path / 'file', content(size))
    monkeypatch.setattr(qh, 'getsize', lambda _: claimed)
    with pytest.raises(OSError, match='Unexpected end of file'):
        qh.quickhash(path, N)
