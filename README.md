# MoDEL Workflow



## Installation

It is recommended to create a virtual environment. You can do this using Conda and the `environment.yml` file, with:

```conda env create --file environment.yml```

From the workflow's directory then execute the following:

```python setup.py install```

This will add the package directory to `PATH` and create the entry point `mwf` to be used when executing the workflow.

---

## How to use

print help message:

```mwf --help```

or  ```mwf -h```

**actions**

* run analysis on current directory, which is asumed to contain topology, trajectory and input files

```mwf```

you can also specify a different directory with the input files using the `-dir` or `--working_dir` options.

* run analysis from already uploaded project

```model_wf -p <project name>```

By default the data files will be downloaded from `https://bioexcel-cv19-dev.bsc.es`. Another URL can be specified with the `-url` option.

* run single analysis from workflow

```model_wf -a <analysis name>```

Analysis name choices: `{rmsd,rmsf,rgyr,pca,pca_contacts,rmsd_per_residue,rmsd_pairwise,distance_per_residue,hydrogen_bonds,energies,pockets}`

_Other options:_

`-in`: input filename for the analysis. Default: `inputs.json`

