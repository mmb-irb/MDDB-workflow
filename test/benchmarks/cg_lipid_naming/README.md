# Coarse grained lipids by bead names

Benchmark of `coarse_grain_lipid_tails` in `mddb_workflow/tools/cg_lipids.py`. The membrane mapping uses it to find
coarse grained lipids whose residue names are not in `LIPIDS_RESIDUE_NAMES`, and it also gives the tail beads of lipids
for the membrane mapping and the lipid analyses. The rule only uses bead names, because MDDB topologies have no bead
types and bonds may be missing.

## Files

- `corpus.json`: 1970 unique (residue name, bead names) cases from the parameter sets of Martini and other coarse grained
  force fields, labelled as `fatty`, `sterol` or `null` (not a lipid). Each case also has its class of molecule, the
  source groups where it is found and the first file where it is found.
- `build_corpus.py`: downloads the parameter sets and builds `corpus.json`.
- `test_cg_lipid_naming.py`: checks the rule against the corpus. The lipids it misses and the non-lipids it accepts must
  be the recorded ones, so the test fails on regressions and on improvements alike. Update the record and this README
  when the rule improves.

## Running it

```bash
pytest test/benchmarks/cg_lipid_naming
```

To print the results by source as a markdown table:

```bash
python test/benchmarks/cg_lipid_naming/test_cg_lipid_naming.py
```

To rebuild the corpus, e.g. to add a parameter set or fix a label:

```bash
python test/benchmarks/cg_lipid_naming/build_corpus.py
```

It downloads about 900 MB into `downloads/` next to the script (ignored by git, or use `--cache DIR`) and takes about
two minutes. Git repositories are pinned to commits in `build_corpus.py`. The Martini 2 lipidome files are listed from the
live cgmartini page, so a rebuild may change if the page changes. The committed corpus is the reference.

## Results

The rule recognises 1263 of the 1305 lipids (96.8%) and accepts 9 of the 665 non-lipids (1.4%), as of 2026-10-26.

| source | lipids | recognised | non-lipids | false positives |
|---|---|---|---|---|
| insane | 170 | 166 | 4 | 0 |
| M2 lipidome | 211 | 203 | 37 | 1 |
| M2 other | 25 | 25 | 137 | 5 |
| M3 ionizable-lipid repo | 306 | 306 | 88 | 3 |
| M3 lipidome | 670 | 670 | 96 | 0 |
| M3 official non-lipids | 0 | 0 | 190 | 0 |
| M3 small-molecule library | 0 | 0 | 91 | 0 |
| M3 sterols | 1 | 1 | 0 | 0 |
| M3 v1 phospholipids | 112 | 112 | 0 | 0 |
| polyply/vermouth | 0 | 0 | 147 | 0 |
| SIRAH | 5 | 0 | 7 | 0 |
| SPICA | 29 | 0 | 60 | 0 |
| TS2CG | 125 | 125 | 0 | 0 |
| **Total** | **1305** | **1263** | **665** | **9** |

A case counts in every source where it is found, so the source rows add up to more than the total.

Missed lipids (they must be found by residue name, `LIPIDS_RESIDUE_NAMES`, or by the `forced_lipids` input):
- SPICA (29) and SIRAH (5) lipids, which name their beads in other ways (`C11`, `BC11`...).
- Hopanoids (4): they have no `ROH` head bead, and a looser rule on `R1`, `R2`... beads would accept many small molecules.
- Lyso-PC with one or two tail beads (`CPC`, `TPC`), `BOG` (two tail beads) and the Martini 2 single residue LPS (`V11`... beads).

False positives:
- 5 two tailed voltage-sensitive dyes, which are membrane-embedded amphiphiles anyway.
- 3 parameterisation fragments of ionizable lipids, which never occur in production systems.
- 1 acyl residue of a multi-residue lipid A molecule.

### History

The rule was designed with this corpus on 2026-10-23:
- **v1 prototype:** `^[CDT]\d{1,2}[A-D]$` tail beads, or `ROH` and sterol rings. It recognised 1189 lipids with 17 false
  positives. It missed every cardiolipin, the C8-C10 lipids, the branched ionizable lipids and the sterols with an `R0`
  head.
- **v2 without the single tail extension:** 1245 lipids recognised.
- **v2 with the single tail extension (current):** 1263 lipids recognised and 9 false positives.

The corpus had 755 non-lipids then, including 90 all atom (OPLS) models of the Martini 3 small-molecule library. They
are no longer included, since they are not coarse grained.

## Labels

- **Lipids:**
  - phospholipids, sphingolipids, ceramides and glycolipids;
  - sterols, including sterol esters and glycosides, and hopanoids;
  - DAG, MAG and TAG, fatty acids, detergents and bolalipids;
  - cationic, ionizable and PEG lipids, cardiolipins and the single residue LPS.
- **Not lipids:**
  - pure hydrocarbons (no polar bead) and free inositol phosphates;
  - separate PEG monomer residues of PEG-lipids;
  - parameterisation fragments of ionizable lipids and residues of multi-residue lipid A molecules;
  - membrane dyes;
  - proteins, including lipidated amino acids, nucleic acids, sugars and polymers;
  - solvents, ions, cofactors and small molecules.

`build_corpus.py` sets the labels by file and section. Not checked, because their parameters were not available: Martini
3 LPS, Martini 3 GM1/GM2, Martini 3 DPC and the CHARMM-GUI Martini Maker libraries.

## Sources

- Martini 3 lipidome: https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters (Pedersen et al., ACS Cent. Sci. 2025)
- Martini 3 sterols: https://github.com/Martini-Force-Field-Initiative/M3-Sterol-Parameters (Borges-Araújo et al., JCTC 2023)
- Martini 3 ionizable lipids, PEG-lipids and sterol derivatives: https://github.com/Martini-Force-Field-Initiative/M3_Ionizable_Lipids (Kjølbye et al., JCTC)
- Martini 3 v1 phospholipids, sugars, small molecules, nucleobases, solvents and ions: https://github.com/marrink-lab/martini-forcefields
- Martini 3 small molecules library: https://github.com/ricalessandri/Martini3-small-molecules
- Martini 2 lipidome and other molecules (polymers, solvents, dry Martini, DNA, RNA, cofactors, dyes): https://cgmartini.nl/docs/downloads/force-field-parameters/martini2/lipidome.html and its sibling pages
- insane lipid templates: https://github.com/Tsjerk/Insane (Wassenaar et al., JCTC 2015)
- TS2CG Martini 3 lipid library and cardiolipins: https://github.com/marrink-lab/TS2CG1.1
- polyply and vermouth (polymers, proteins, DNA): https://github.com/marrink-lab/polyply_1.0 and https://github.com/marrink-lab/vermouth-martinize
- SPICA: https://github.com/SPICA-group/spica-tools (Seo & Shinoda, JCTC 2019)
- SIRAH: https://github.com/SIRAHFF/documentation (GROMACS release `sirah_x2.2_20-07.ff`; Barrera, Machado, Pantano, JCTC 2019)
