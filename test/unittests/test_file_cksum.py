"""Tests of the file cksum: the file size and its quickhash (the whole file, or 3 samples of the large files if unsafe)."""
import pytest
import xxhash

from mddb_workflow.tools.viewer import FILE_CKSUM_PATTERN
from mddb_workflow.utils import file as file_module
from mddb_workflow.utils.file import File
from mddb_workflow.utils.quickhash import quickhash

pytestmark = pytest.mark.CI


def write(path, data: bytes) -> File:
    path.write_bytes(data)
    return File(str(path))


def test_small_files_keep_their_cksum(tmp_path):
    """Files smaller than 3 samples keep the xxh64 of their whole content, so the caches remain valid
    (2.5 MiB: read in 3 chunks, the same value as reading the file at once)."""
    for data in (b'', bytes(range(256)) * 10240):
        for unsafe in (False, True):
            assert write(tmp_path / 'small.xtc', data).get_cksum(unsafe=unsafe) == f'{len(data)}-{xxhash.xxh64(data).hexdigest()}'


def test_large_files(tmp_path, monkeypatch):
    """Files of at least 3 samples (16 bytes here): sampled if unsafe, otherwise hashed whole."""
    monkeypatch.setattr(file_module, 'QUICKHASH_SAMPLE_SIZE', 16)
    data = bytes(range(256)) * 20
    path = tmp_path / 'large.xtc'
    cksum = write(path, data).get_cksum(unsafe=True)
    assert cksum == f'{len(data)}-(UNSAFE){quickhash(str(path), 16)}'
    assert cksum != f'{len(data)}-(UNSAFE){xxhash.xxh64(data).hexdigest()}'
    # The size is still read by the viewer to find the cached files
    assert FILE_CKSUM_PATTERN.match(cksum).group(1) == str(len(data))
    # The safe cksum reads the whole file
    assert write(path, data).get_cksum() == f'{len(data)}-{xxhash.xxh64(data).hexdigest()}'


def legacy_cksum(data: bytes, limit: int = None) -> str:
    """The cksum as it was calculated before using the quickhash."""
    if limit is not None and len(data) > limit:
        return f'{len(data)}-(UNSAFE){xxhash.xxh64(data[:limit]).hexdigest()}'
    return f'{len(data)}-{xxhash.xxh64(data).hexdigest()}'


def test_legacy_cksum(tmp_path, monkeypatch):
    """The legacy cksum is the xxhash of the whole file, or of its head if unsafe and large."""
    monkeypatch.setattr(File, 'LEGACY_CKSUM_UNSAFE_SIZE_LIMIT', 100)
    data = bytes(range(256)) * 20
    path = tmp_path / 'large.xtc'
    assert write(path, data).get_legacy_cksum(unsafe=True) == legacy_cksum(data, 100)
    assert write(path, data).get_legacy_cksum() == legacy_cksum(data)
    assert write(path, data[:50]).get_legacy_cksum(unsafe=True) == legacy_cksum(data[:50], 100)


def test_legacy_cksums_in_cache_match(tmp_path, monkeypatch):
    """Cached legacy cksums are accepted, also inside lists of files, but changed files are not."""
    from mddb_workflow.utils.arg_cksum import get_cksum_id, cksum_matches
    monkeypatch.setattr(file_module, 'QUICKHASH_SAMPLE_SIZE', 16)
    monkeypatch.setattr(File, 'LEGACY_CKSUM_UNSAFE_SIZE_LIMIT', 100)
    data = bytes(range(256)) * 20
    path = tmp_path / 'large.xtc'
    value = [write(path, data), 'other']
    new_cksum = get_cksum_id(value)
    old_cksum = get_cksum_id([write(path, data), 'other'], legacy=True)
    assert new_cksum != old_cksum
    assert cksum_matches(value, new_cksum, old_cksum)
    assert cksum_matches(value, new_cksum, new_cksum)
    assert not cksum_matches(value, new_cksum, None)
    # A change in the file makes both cksums miss
    changed = [write(path, b'\xff' + data[1:]), 'other']
    assert not cksum_matches(changed, get_cksum_id(changed), old_cksum)
