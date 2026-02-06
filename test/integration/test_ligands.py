from mddb_workflow.utils.file import File
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.structures import Structure
from mddb_workflow.tools.get_inchi_keys import InChIKeyData, generate_inchi_references
from mddb_workflow.tools.get_ligands import generate_ligand_references
from mddb_workflow.tools.residue_mapping import generate_residue_mapping
import pytest
import pathlib


@pytest.mark.unit_int
def test_generate_ligand_references():
    """Test the generate_ligand_references function."""
    pdb_file = pathlib.Path(__file__).parent.parent / 'data/input/structures/cin_A000V_structure.pdb'

    mwf_stc = Structure.from_file(str(pdb_file))
    cache = Cache(File('cache_ligands.json'))
    # TODO: make test of forced selection of a ligands from the inputs.yaml
    input_ligands = [
        '5957',
        {'drugbank': 'DB00945'},
        {'chembl': 'CHEMBL112', 'name': 'Paracetamol', 'selection': 'resname LIG'}]  # 'RZVAJINKPMORJF-UHFFFAOYSA-N'
    pdb_ids = ['2L94']
    test_cases = [
        # Direct InChIKey matching
        ('MPUVBVXDFRDIPT-CHKWXVPMSA-N',
         'InChI=1S/C8H13NO2/c9-8(7(10)11)4-5-1-2-6(8)3-5/h5-6H,1-4,9H2,(H,10,11)/t5-,6+,8+/m1/s1',
         ['13885838']),
        # After neutralization
        ('ICACWKLCNCPHAM-YJBOKZPZSA-L',
         'InChI=1S/C21H25N2O8P/c1-15(23-21(27)30-13-17-10-6-3-7-11-17)19(24)22-14-32(28,29)31-18(20(25)26)12-16-8-4-2-5-9-16/h2-11,15,18H,12-14H2,1H3,(H,22,24)(H,23,27)(H,25,26)(H,28,29)/p-2/t15-,18-/m0/s1',
         ['16741275']),
        # After removing stereochemistry
        ('GYDJEQRTZSCIOI-LJGSYFOKSA-N',
         'InChI=1S/C8H15NO2/c9-5-6-1-3-7(4-2-6)8(10)11/h6-7H,1-5,9H2,(H,10,11)/t6-,7-',
         ['5526']),
        # After neutralization + stereochemistry removal
        ('CLRSLRWKONPSRQ-IIPSPAQQSA-O',
         'InChI=1S/C38H47ClN4O4/c1-25(2)47-35-22-33-28(20-34(35)46-5)21-36(44)43(38(33)27-8-10-29(39)11-9-27)32-16-14-30(15-17-32)41(4)23-26-6-12-31(13-7-26)42-19-18-40(3)37(45)24-42/h8-11,14-17,20,22,25-26,31,38H,6-7,12-13,18-19,21,23-24H2,1-5H3/p+1/t26-,31-,38-/m0/s1',
         ['58437678']),
        # After standardization
        ('ZGIPMGGESRKDMZ-XRNKAMNCSA-N',
         'InChI=1S/C22H23ClN2O8/c1-21(32)7-6-8-15(25(2)3)17(28)13(20(24)31)19(30)22(8,33)18(29)11(7)16(27)12-10(26)5-4-9(23)14(12)21/h4-5,7-8,15,26,29-30,32-33H,6H2,1-3H3,(H2,24,31)/t7-,8-,15-,21-,22-/m0/s1',
         ['54675777']),
        # PDB fallback
        ('UZESPYYUUQTAQR-IQRFGFHNSA-P',
         'InChI=1S/C18H34N6/c1-23(2)13-5-11-21-17(19)15-7-9-16(10-8-15)18(20)22-12-6-14-24(3)4/h7-10,21-22H,5-6,11-14,19-20H2,1-4H3/p+2/b17-15-,18-16-',
         ['53245673', '137349678']),
        # From DrugBank
        ('BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
         'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
         ['2244']),
        # From inputs.yaml
        ('ZKHQWZAMYRWXGA-KQYNXXCUSA-J',
         'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1',
         ['5957']),  # Workflow finds 5461108, Atp(4-)

    ]
    inchimap = {case[0]: InChIKeyData(inchi=case[1], resnames=set([1])) for case in test_cases}
    lipid_references = {}
    memmap = {'no_mem_lipid': []}

    ligand_references = generate_ligand_references(
        mwf_stc, cache,
        input_ligands,
        pdb_ids,
        inchimap,
        lipid_references,
        memmap
    )
    assert 'RZVAJINKPMORJF-UHFFFAOYSA-N' in ligand_references, "Forced CHEMBL112 not found"
    assert len(ligand_references) == len(test_cases) + 1
    for case in test_cases:
        found = ligand_references[case[0]]['pubchem']
        assert found in case[2], f"Ligand {case[0]} found with CID {found} instead of {case[2]}."

    inchirefs = generate_inchi_references(inchimap, {}, ligand_references, File('inchiout.json'))
    residue_map = generate_residue_mapping([], inchirefs, mwf_stc)
    assert 'ZKHQWZAMYRWXGA-KQYNXXCUSA-N' in residue_map['references'],\
        'ZKHQWZAMYRWXGA-KQYNXXCUSA-J should be transformed to ZKHQWZAMYRWXGA-KQYNXXCUSA-N'
    # Cleanup
    pathlib.Path('cache_ligands.json').unlink()
    pathlib.Path('inchiout.json').unlink()


if __name__ == '__main__':
    test_generate_ligand_references()
