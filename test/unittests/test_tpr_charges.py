"""Tests for the TPR charges mined with Gromacs, used when MDAnalysis can not read the TPR (e.g. tpx 102, FAIRMD 472)."""
import os
import pytest

from mddb_workflow.utils.gmx_spells import get_tpr_atom_count, get_tpr_charges

# The dummy Ala-Ala-Ala peptide twice, 2 TIP3P waters, Na+ and Cl- (74 atoms in 4 molecule types)
REPEATED_MOLECULES_TPR = 'input/dummy/gromacs/repeated_molecules.tpr'


@pytest.mark.unit_int
def test_tpr_charges_cover_every_atom(test_data_dir: str, tmp_path, monkeypatch):
    """Every atom gets its charge, also in molecule types with several molecules."""
    tpr_filepath = os.path.join(test_data_dir, REPEATED_MOLECULES_TPR)
    monkeypatch.chdir(tmp_path)
    charges = get_tpr_charges(tpr_filepath)
    assert len(charges) == get_tpr_atom_count(tpr_filepath) == 74
    # Both peptides have the same charges
    assert charges[:33] == charges[33:66]
    # TIP3P water, Na+ and Cl-
    assert charges[66:] == pytest.approx([-0.834, 0.417, 0.417] * 2 + [1, -1])
    assert sum(charges) == pytest.approx(0, abs=1e-4)
