# Generating inputs.yaml files

Each project needs an `inputs.yaml` before the workflow runs. `mwf dataset inputs` (or
`Dataset.generate_inputs_yaml`) renders one Jinja2 **template** per project, optionally fed by a
**generator** that computes per-project values.

> **Copy the bundled example first.** `assets/inputs_template.yaml` and `assets/nanobodies_generator.py`
> are a complete, working template/generator pair (a GROMACS nanobody dataset — the exact layout in the
> worked example below). For a new dataset, copy them and edit the metadata rather than starting from
> scratch. The canonical, fully-commented reference template lives at
> `mddb_workflow/resources/inputs_file_template.yml`.

## The one rule that must not be skipped

**Every template must set `dataset_path: {{DATASET}}`.** `DATASET` renders to the dataset file's path
relative to each project, so the generated `inputs.yaml` points back to the dataset the project belongs
to and the workflow updates the correct dataset's state. Omitting it, or hardcoding an absolute path,
breaks status tracking. Put it near the top of the template and treat it as non-negotiable:

```yaml
dataset_path: {{DATASET}}
```

## First understand the folder structure — it drives the input file fields

Before writing a template or generator, **inspect a representative project directory and map its files
to the four input-file fields.** These fields tell the workflow what to actually read, so getting them
right matters more than any metadata — a wrong path here means the workflow reads the wrong data or
fails. Never assume a naming convention; list the real files (`ls -R` a sample project) and confirm.

| Field | What it is | Cardinality |
|-------|-----------|-------------|
| `input_topology_filepath` | The topology (bonds/charges) — a GROMACS `.tpr`/`.top`, Amber `.prmtop`, etc. | **Exactly one per project.** A structure may substitute only if topology is missing (set `input_topology_filepath: no` and provide a structure). |
| `input_structure_filepath` | Starting coordinates (`.pdb`/`.gro`). Optional when a self-contained topology like a `.tpr` is given; **required** when topology is absent. Must match the trajectory's atom set (e.g. the *solvated* structure if the trajectory is solvated). | At most one per project. |
| `input_trajectory_filepaths` | The trajectory file(s). Can be set project-wide, but **it almost always belongs per-MD**, since each replica has its own. | Per-MD (preferred) or project-level. |
| `mds` | The list of MDs (replicas/subprojects). Each entry needs an `mdir` (the folder that MD lives in) and typically its own `input_trajectory_filepaths`. Per-MD keys override project metadata. | One list; one entry per replica. |

**`mdir` must match the real directory name on disk.** This is the most common mistake: writing
`mdir: replica_1` when the folder is actually called `r1`. The workflow uses `mdir` to locate the MD's
files, so a mismatch silently points it at the wrong (or missing) directory. This is the single biggest
reason to derive `mds` from the filesystem in the generator rather than hand-writing them.

### Worked example: a GROMACS nanobody project

```
2bse_F2-123_NA_run/                     ← project directory
├── 2bse_F2-123_NA_solv_ions.gro        ← solvated+ionised structure (matches the trajectory)
├── topol.top                           ← project topology
├── em.* nvt.* npt.* ions.*             ← preparation/equilibration (NOT production inputs)
├── r1/  md.xtc  md.tpr  md.gro  topol.top   ← replica 1: production trajectory + topology
├── r2/  md.xtc  md.tpr  ...                 ← replica 2
└── r3/  md.xtc  md.tpr  ...                 ← replica 3
```

The correct mapping, and the reasoning:

- **Topology** → a `.tpr` is self-contained (topology + coordinates), so the most robust choice is one
  representative production `r*/md.tpr`; the replicas share the same molecular topology. (Using the
  project-level `topol.top` also works but relies on its `.itp` includes being resolvable.)
- **Structure** → the solvated `*_solv_ions.gro`, because it matches the atom set in the solvated
  `md.xtc`. The bare `.pdb` is protein-only and would *not* match — a concrete example of why you must
  check, not guess. Optional if a `.tpr` topology is used.
- **Trajectories** → per-MD: `r1/md.xtc`, `r2/md.xtc`, `r3/md.xtc`. The many prep files (`em`, `nvt`,
  `npt`) are equilibration, not production — do not feed them as trajectories.
- **`mds`** → three entries with `mdir: r1`, `r2`, `r3` (the real folder names), each carrying its own
  `input_trajectory_filepaths`.

The generator's job is to discover exactly this from disk. `assets/nanobodies_generator.py` and
`assets/inputs_template.yaml` are that generator/template pair for this very layout — the point is that
the generator encodes *your* directory structure, so adapt the globs when a new dataset differs.

## The template

Start from the canonical, fully-commented template at
`mddb_workflow/resources/inputs_file_template.yml` and trim it to what the dataset needs. The template
is rendered once per project and always receives these variables:

- `DIR` — the project directory name (e.g. `project_1`).
- `PATH` — the full project path (a `Path`).
- `DATASET` — the dataset file path, relative to the project directory.
- `project_status` — a dict of the project's DB row (`uuid`, `rel_path`, `num_mds`, `state`, …).
- Everything the generator returns (merged in as top-level variables).

Example template:

```yaml
name: Project {{DIR}}
dataset_path: {{DATASET}}
description: 10 ns simulation of {{DIR}}
{%- if is_special_case %}
authors:
  - Special Case Author
{%- else %}
authors:
  - Default Author
{%- endif %}

mds:
{% for md in mds %}
  - mdir: {{ md.mdir }}
    input_trajectory_filepaths: {{ md.traj }}
{% endfor %}
```

Note `is_special_case`, `mds`, `md.mdir`, `md.traj` are **not** built in — they come from the
generator below. The template and generator are a matched pair: the template consumes exactly the keys
the generator produces.

## The generator

A `.py` file that defines a function with this exact signature:

```python
def inputs_generator(project_dir: str) -> dict | None:
    ...
```

It is called once per project with the project directory path and returns a dict of variables for the
template (or `None`/empty for defaults). Keep it pure and fast — it runs for every matched project.

Simple flag example:

```python
def inputs_generator(project_dir: str):
    if "special_cases" in project_dir:
        return {'is_special_case': True}
```

Discovering MDs from files on disk:

```python
from pathlib import Path

def inputs_generator(project_dir: str):
    """Build the MD list from production trajectory files in the project."""
    project_path = Path(project_dir)
    prod_trajs = sorted(project_path.glob('prod_*.traj'))
    mds = [
        {
            'mdir': f'md_replica_{i + 1}',
            'traj': traj.relative_to(project_path).as_posix(),
        }
        for i, traj in enumerate(prod_trajs)
    ]
    return {'mds': mds}
```

## Running it

```bash
# Template only (uses template defaults for everything)
mwf dataset inputs mwf_ds.sql -it template.yaml

# Template + generator, restricted to new projects, overwriting existing files
mwf dataset inputs mwf_ds.sql -it template.yaml -ig generator.py -o -st new
```

Without `-o/--overwrite`, projects that already have an `inputs.yaml` are skipped and reported — so
re-running is safe and only fills in the missing ones. Pass `-o` only when the user wants to replace
existing files.

Programmatic equivalent (pass a callable directly instead of a path):

```python
ds.generate_inputs_yaml(
    inputs_template_path='template.yaml',
    inputs_generator=inputs_generator,   # a callable, or a path string to the .py file
    overwrite=False,
    query_state=['new'],
)
```

## Pre-flight: dry-run the generator over the whole dataset first

Before committing a large dataset to the cluster, **run the generator against every project and
inspect the results** — it's fast, read-only, needs no `mwf`, and catches problems that would otherwise
surface only as failed jobs (or, worse, as silently wrong inputs). A generator maps the filesystem to
input fields, so it's the perfect place to spot layout anomalies: projects still copying, missing
trajectories, an unexpected number of replicas, or a structure/topology file that isn't where you
assumed. This check has real teeth — on the nanobody dataset it revealed that 16 of 220 projects had
two replicas rather than three, and one had an `r3/` folder with no `md.xtc` (an incomplete replica the
generator correctly dropped). None of that is visible from a single spot-checked project.

Two things to look at: (1) does the generator resolve for every project, and (2) is the distribution of
what it found what you expected (e.g. replica counts). A silent, uniform result is reassuring; an
outlier is a project to inspect before launching.

```python
from pathlib import Path
from collections import Counter
from my_generator import inputs_generator   # import your generator

projects = sorted(p for p in Path('vhh_run').iterdir() if p.is_dir())
dist, failed = Counter(), []
for p in projects:
    try:
        g = inputs_generator(str(p))
        dist[len(g['mds'])] += 1                     # e.g. replica-count distribution
        # sanity-check that the discovered paths actually exist
        for md in g['mds']:
            assert (p / md['traj']).exists(), f"missing trajectory in {p.name}"
    except Exception as e:
        failed.append((p.name, str(e)))

print(f"{len(projects)} projects | replica counts: {dict(sorted(dist.items()))}")
print("FAILED:", failed or "none")
```

Investigate any failures and any minority buckets in the distribution. Decide what to do with the
oddballs (skip an abandoned replica, wait for a copy to finish, fix the layout) *before* you generate
inputs and launch — that is far cheaper than debugging it from job logs afterward. Note that a generator
that raises on incomplete projects (rather than silently returning empty `mds`) makes this check far
more useful, so prefer raising over returning nothing when required files are absent.
