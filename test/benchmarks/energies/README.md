# Interaction energies benchmark

Benchmark of the energies analysis (`mddb_workflow/analyses/energies.py`): the electrostatic and van der Waals
interaction energy of every atom of an agent with the other agent of an interaction. It compares the analysis' two
engines, chosen with `mwf run --energies_engine`:

- `cmip` (the default): CMIP, which solves the Poisson-Boltzmann equation on a 0.5 Å grid.
- `gridfree`: Coulomb with the Mehler-Solmajer sigmoidal distance-dependent dielectric (bulk 78.4, no cutoff) and
  CMIP's van der Waals table (radii scaled by 0.8, 8 Å cutoff), summed over atom pairs in a numba kernel.

Every system goes through the workflow's own `energies()` function, so the benchmark measures the code users run. An
engine is good if its per-atom energies rank experimental binding data: the energy of a ligand against the measured
affinities of a congeneric series, and the energy of a mutated residue's side chain against the measured effect of the
mutation to alanine.

## Results

Scores of both engines on the same inputs (energies analysis 1.2.0, runs of 2026-11-13; SKEMPI-40 rescored on
2026-11-16, when `skempi_prepare.py` kept the sign of negative residue numbers and 9 mutations of 4P23 joined the
score). Mean Spearman correlation per target, series or complex; the difference is grid-free minus CMIP, paired by
target, series or complex, with its 95% bootstrap interval:

| Benchmark | Units | CMIP | Grid-free | Difference [95% CI] |
|---|---|---|---|---|
| FEP+ (T3) | 37 series, 725 complexes | 0.424 | 0.463 | +0.039 [−0.016, +0.091] |
| PL-REX (T2) | 10 targets, 164 complexes | 0.386 | 0.459 | +0.073 [+0.010, +0.150] |
| SKEMPI-40 (T3) | 40 complexes, 1138 mutations | 0.407 | 0.456 | +0.049 [+0.018, +0.081] |
| SKEMPI-40 hot-spot AUC | | 0.763 | 0.772 | +0.009 [−0.012, +0.030] |
| SKEMPI-7 (T2 part of SKEMPI-40) | 7 complexes, 376 mutations | 0.431 | 0.413 | −0.018 [−0.070, +0.027] |
| Protein-nucleic acid (T4) | 39 complexes, 355 mutations | 0.202 | 0.216 | +0.014 [−0.028, +0.055] |
| Protein-nucleic acid hot-spot AUC | 17 complexes | 0.741 | 0.728 | −0.013 [−0.036, +0.009] |

Time of the analysis, summed over the systems of a tier (one frame each, a busy 16-core laptop):

| Tier | Systems | CMIP | Grid-free | Median per system |
|---|---|---|---|---|
| T3 (includes T2) | 929 | 11.1 h | 4.6 min | 38 s against 0.28 s |
| T4 | 39 | 42 min | 13 s | 49 s against 0.30 s |

For reference, the explicit-water FEP+ calculations of the FEP+ benchmark rank the same 37 series at 0.66, and CMIP
1.1.0 (before the fixes of 1.2.0) scored 0.368 on PL-REX and 0.395 on SKEMPI-7.

T1, one frame each with one OpenMP thread: total energies of each agent in kcal/mol (electrostatic, then van der
Waals), as `t1` prints them. With grid-free both directions of each interaction give the same electrostatic total;
with CMIP they are up to 0.6% apart.

| System | Engine | es agent 1 | es agent 2 | vdw agent 1 | vdw agent 2 | Time (s) |
|---|---|---|---|---|---|---|
| A01TW | CMIP | -7.90 | -7.90 | -20.00 | -20.07 | 55 |
| A01TW | grid-free | -4.33 | -4.33 | -20.72 | -20.72 | 0.5 |
| A0221 | CMIP | -49.72 | -49.73 | -15.89 | -15.90 | 136 |
| A0221 | grid-free | -26.49 | -26.49 | -16.19 | -16.19 | 0.5 |
| 1VFB | CMIP | -30.32 | -30.14 | -26.85 | -26.65 | 67 |
| 1VFB | grid-free | -9.61 | -9.61 | -27.74 | -27.74 | 0.3 |

## Tiers

- **T0**, synthetic: `test/unittests/test_energies.py`, a 37-atom system with golden outputs for both engines.
- **T1**, regression: one frame of three real systems, compared per atom with an earlier run.
  - `A01TW`: kinase-crizotinib, 4850 atoms (MDDB, frame 3 of a 5-frame sample).
  - `A0221`: TNF-antibody, 8915 atoms (MDDB, frame 2 of a 5-frame sample).
  - `1VFB`: SKEMPI D1.3 Fv-lysozyme, 5364 atoms (crystal structure, PDBFixer, Amber ff14SB charges).
- **T2**, small accuracy: PL-REX and the 7 original SKEMPI complexes.
- **T3**, full accuracy: T2 plus the FEP+ benchmark and 33 more SKEMPI complexes. A T3 run computes only the extra
  systems; its scores also use the T2 results of the same label.
- **T4**, protein-nucleic acid alanine scanning (mmCSM-NA), scored like SKEMPI.
- **T5**, sealed hold-out: the 49 SKEMPI alanine-scan complexes not in T2 or T3 (853 mutations; 844 before the
  residue-number fix of 2026-11-16). It is kept for final decisions: `t5 score` needs `--decision` and records each
  use in `downloads/t5/ledger.jsonl`. It has been scored once (2026-11-13, a desolvation term added to grid-free, not
  adopted); the grid-free engine itself has not been scored on it.

Scores: for PL-REX and FEP+, the Spearman correlation between the ligand's total energy and the experimental dG of each
target or series (5 complexes or more). For SKEMPI and T4, the Spearman correlation between the negated side-chain
energy of each mutated residue (both agents) and the experimental ddG of its mutation to alanine, per complex, and the
AUC of hot spots (ddG ≥ 2 kcal/mol). Each target, series or complex counts once in the mean.

## Files

- `__main__.py`: the entry point, `python -m test.benchmarks.energies <tier> ...`.
- `t1.py`: the T1 runner, and the function that runs one frame through `energies()` for every tier.
- `t2.py`: the T2 to T5 runner and scorer.
- `t1_inputs.py`: fetches and builds the T1 frames (coordinates, atoms, residues, chains, charges, bonds and agents).
- `scripts/`: the scripts that prepare the inputs of each tier (see below).
- `downloads/` (ignored by git): sources, prepared inputs, per-atom results (`downloads/<tier>/results/<label>/`) and
  scores (`downloads/<tier>/scores/<label>.json`).

## Replicating

### Environments

- `mwf_env` (`envs/environment.yml`): the workflow, CMIP 2.7.0, numba 0.68.1, AmberTools 24.8 (antechamber, tleap)
  and RDKit 2024.03.6. Runs every tier, `plrex_ligands.sh`, `plrex_build.sh` and `fep_ligand.sh`.
- An OpenMM environment for the structure preparation:
  `mamba create -n openmm_env -c conda-forge python=3.11 openmm=8.6.1 pdbfixer=1.12 numpy scipy`.
  `fep_minimize.py` needs CUDA; the other scripts run on the CPU (`OPENMM_DEFAULT_PLATFORM=CPU`).

### Inputs

| Data | Source | Version |
|---|---|---|
| PL-REX | `https://github.com/Honza-R/PL-REX` | commit `c5a84f8d481a4c811ea2da51c9395b75070acea9` |
| FEP+ benchmark | `https://github.com/schrodinger/public_binding_free_energy_benchmark` | commit `a65d86b782a678c75388566ce00446e45283ee63` |
| SKEMPI 2.0 | `https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv` | sha256 `76f8e60d09eaa4f05fb265063886859a86895c0b2fb225e31fe82928a883972d` |
| SKEMPI and T4 structures | `https://files.rcsb.org/download/<PDB>.pdb` (downloaded by the scripts) | not versioned (2026-11-12 and 2026-11-13) |
| mmCSM-NA single mutations | `https://biosig.lab.uq.edu.au/mmcsm_na/static/datasets/single_mutations.csv` | sha256 `f829c4200a2c8395ae8ca0bbd7be631be12752668a2f151d46942e019a2bba08` |
| T1 frames | MDDB A01TW and A0221 (irb-dev node, the workflow's default API), SKEMPI 1VFB_AB_C | not versioned (fetched frames equal those of 2026-11-12 on 2026-11-16) |

From the repository root, with the sources downloaded to `$SOURCES` (absolute paths: the shell scripts change
directory):

```bash
B=$PWD/test/benchmarks/energies; S=$B/scripts; D=$B/downloads
```

**T2, PL-REX** (`downloads/t2/plrex/<target>/<complex>/atoms.npz` and each target's `experimental_dG.txt`). In
`mwf_env`, ligand parameters (GAFF2, AM1-BCC) and Amber topologies of each complex; then, in the OpenMM environment,
`plrex_mm.py` writes `atoms.npz` (and MM energies, not used here):

```bash
for c in $SOURCES/PL-REX/0*/structures_pl-rex/*/; do t=$(basename $(dirname $(dirname $c))); n=$(basename $c); $S/plrex_ligands.sh $c/ligand.sdf $D/t2/plrex/$t/$n && $S/plrex_build.sh $c $D/t2/plrex/$t/$n; done
```

```bash
for t in $SOURCES/PL-REX/0*/; do cp $t/experimental_dG.txt $D/t2/plrex/$(basename $t)/; done
```

```bash
python $S/plrex_mm.py $D/t2/plrex $D/t2/plrex_mm.json
```

**T2 and T3, SKEMPI** (`downloads/<tier>/skempi/<complex>.npz` and `_ddg.json`), in the OpenMM environment:

```bash
python $S/skempi_prepare.py $SOURCES/skempi_v2.csv $D/t2/skempi 1A22_A_B 1EMV_A_B 1JRH_LH_I 1JTD_A_B 1JTG_A_B 1VFB_AB_C 3S9D_A_B
```

```bash
python $S/skempi_prepare.py $SOURCES/skempi_v2.csv $D/t3/skempi 1A4Y_A_B 1AO7_ABC_DE 1BP3_A_B 1BRS_A_D 1DAN_HL_UT 1DQJ_AB_C 1DVF_AB_CD 1F47_A_B 1K8R_A_B 1KTZ_A_B 1LFD_A_B 1OGA_ABC_DE 1Z7X_W_X 2G2U_A_B 2NYY_DC_A 2WPT_A_B 3BN9_B_CD 3C60_CD_AB 3EQS_A_B 3EQY_A_C 3HFM_HL_Y 3M62_A_B 3M63_A_B 3NGB_HL_G 3QHY_A_B 3QIB_ABP_CD 3SZK_AB_C 4G0N_A_B 4OFY_A_D 4P23_CD_AB 4RS1_A_B 5C6T_HL_A 5M2O_A_B
```

**T3, FEP+** (`downloads/t3/fep/<series>/<ligand>/atoms.npz` and each series' `series.json`). `fep_select.py`
chooses the series and ligands by rules fixed before any energy was computed (10 ligands or more with a structure, an
experimental range of 2 kcal/mol or more, at most 25 ligands per series, seed 0); `fep_ligand.sh` (in `mwf_env`)
parameterizes each ligand; `fep_receptor.py` builds the receptor and the complexes; `fep_minimize.py` relaxes each
complex with its heavy atoms restrained to the input pose, since the FEP+ poses are unrelaxed and many clash:

```bash
python $S/fep_select.py $SOURCES/public_binding_free_energy_benchmark $D/t3/fep
```

```bash
for l in $D/t3/fep/*/lig*/; do $S/fep_ligand.sh $l; done
```

```bash
python $S/fep_receptor.py $D/t3/fep && python $S/fep_minimize.py $D/t3/fep
```

**T4** (`downloads/t4/nucleic/`), in the OpenMM environment:

```bash
python $S/nucleic_prepare.py $SOURCES/single_mutations.csv $D/t4/nucleic
```

**T5** (`downloads/t5/skempi/`): `select` writes the candidates by the tier's rule, `skempi_prepare.py` prepares them
(OpenMM environment) and `check` keeps those with 6 or more mutated residues in the structure. It should select 49
complexes with 853 mutations:

```bash
python $S/t5_select.py select $SOURCES/skempi_v2.csv $D/t5/skempi
```

```bash
OPENMM_DEFAULT_PLATFORM=CPU python $S/skempi_prepare.py $SOURCES/skempi_v2.csv $D/t5/skempi <the candidates printed by select>
```

```bash
python $S/t5_select.py check $D/t5/skempi
```

**T1** (`downloads/t1/<system>.npz`): the runner builds each frame the first time it needs it (`t1_inputs.py`).
For A01TW and A0221 it downloads a 5-frame sample with the workflow (`mwf run -proj <accession> -i download -smp 5`)
and runs its interactions and charges tasks (`mwf run -mdir <replica> -i inter charges -t`), about a minute per
project. 1VFB is the T2 SKEMPI complex 1VFB_AB_C, so prepare the T2 SKEMPI inputs first.

### Running

In `mwf_env`, from the repository root. CMIP runs need several GB; run them in their own memory-capped unit, for
example `systemd-run --user --collect --unit=energies-t3 -p MemoryMax=10G -p MemorySwapMax=0 --working-directory=$PWD
bash -c "..."`.

T1 prints the totals of the table above (CMIP about 4 minutes with one thread, grid-free about 2 s, after the first
fetch). To check that a change keeps the numbers, store a run of the default engine (CMIP) under a label before the
change and compare with it after; the comparison fails if any per-atom value differs by more than 1e-6 kcal/mol:

```bash
python -m test.benchmarks.energies t1 --label before
```

```bash
python -m test.benchmarks.energies t1 --baseline before
```

`--engine gridfree` does the same with the grid-free engine. `--repeat 5 --threads 16` runs each system 5 times and
reports how much the results vary.

T2 to T4 run every system of a tier under a label and score it, with paired differences against another label:

```bash
python -m test.benchmarks.energies t2 run gridfree-1.2.0 --engine gridfree
```

```bash
python -m test.benchmarks.energies t2 run cmip-1.2.0 --engine cmip
```

```bash
python -m test.benchmarks.energies t2 score gridfree-1.2.0 --baseline cmip-1.2.0
```

and the same with `t3` (after `t2`, whose results it reuses) and `t4`. `run` skips the systems a label already has,
so an interrupted run resumes. `score` prints the metrics, to compare with the tables above, and writes them to
`downloads/<tier>/scores/<label>.json`.

T5 is scored only to settle a decision whose rule is fixed beforehand:

```bash
python -m test.benchmarks.energies t5 score <label> --baseline <label> --decision '<the decision and its rule>'
```

### Notes

- CMIP needs about 220 bytes per grid point: five SKEMPI complexes of T3 (2NYY, 3M62, 3M63, 4OFY, 5C6T) need 5-7 GB
  each. Run CMIP on T3 with one job (`--jobs 1 --threads 8`) under a 10 GB cap, or with more jobs under a 5 GB cap and
  `-p OOMPolicy=continue`, and rerun the failures with one job.
- CMIP truncates file paths longer than 150 characters, so the runner works in a short temporary directory.
- On a busy machine CMIP is faster with 8 threads than with 16. Grid-free results do not depend on the thread count.
- Preparing a complex again does not give the same atoms PDBFixer adds: 4P23, prepared again on 2026-11-16 with the
  same versions, had all its hydrogens and 181 heavy atoms of incomplete side chains elsewhere (the other heavy atoms
  were identical). A new preparation therefore gives slightly different energies and scores than the tables above.
- A CMIP rerun of 4P23 with 2 threads differed from the stored run by 0.001 kcal/mol, the last digit CMIP prints,
  in 4% of the atoms. T1 runs CMIP with one thread so that its results are reproducible.
