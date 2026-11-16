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
