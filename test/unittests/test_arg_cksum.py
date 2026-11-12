from mddb_workflow.utils.arg_cksum import get_cksum_id
from mddb_workflow.utils.selections import Selection


def test_selection_cksum_collision():
    """Selections with the same length and index sum must have different cksums."""
    a = Selection(list(range(0, 10)) + list(range(30, 40)))
    b = Selection(list(range(10, 30)))
    assert get_cksum_id(a) != get_cksum_id(b)


def test_selection_cksum_reproducible():
    """The same atom indices give the same cksum, whatever their input order."""
    assert get_cksum_id(Selection([3, 1, 2])) == get_cksum_id(Selection([1, 2, 3]))


def test_set_cksum():
    """Sets are cksummed by their sorted values, not all as None."""
    assert get_cksum_id({1, 2}) is not None
    assert get_cksum_id({1, 2}) != get_cksum_id({1, 3})
    assert get_cksum_id({2, 1}) == get_cksum_id({1, 2})
