"""Tests for the Cache class."""

from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.file import File


def test_parallel_runs_keep_each_other_changes(tmp_path):
    """Test that two runs using the same cache at the same time do not erase each other changes.

    Bug scenario: Every run loads the cache in memory and writes all of it every time something changes.
    Before writing, the cache file should be read again if another run changed it.
    However, the cksum of the cache file was calculated only once, so changes from other runs were never seen.
    Then the last run writing the cache erased the changes of the other run.
    """
    cache_filepath = str(tmp_path / '.mwf_cache.json')
    first_run_cache = Cache(File(cache_filepath))
    second_run_cache = Cache(File(cache_filepath))

    first_run_cache.update('first', 1)
    second_run_cache.update('second', 2)
    first_run_cache.update('third', 3)

    final_cache = Cache(File(cache_filepath))
    assert final_cache.retrieve('first') == 1
    assert final_cache.retrieve('second') == 2
    assert final_cache.retrieve('third') == 3
    assert first_run_cache.retrieve('second') == 2
