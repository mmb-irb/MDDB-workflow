# MoDEL Workflow



## Installation

This will install the package and create the entry point `mwf` to be used when executing the workflow.

It is recommended to create a virtual environment. You can do this using Conda and the `environment.yml` file, with:

`conda env create --file environment.yml`

From the workflow's directory then execute the following:

In DEVELOPMENT:

`python setup.py develop`

To uninstall the develop mode use:

`python setup.py develop -u`

or:

`pip uninstall model-workflow`

To make notebooks fully operative use:
WARNING: This may break the environment for using the mwf

`conda install -c conda-forge notebook`

In PRODUCTION:

`python setup.py install`

To uninstall the production mode use:

`python setup.py install --record filelist`

`xargs rm -rf < filelist`

---

## Installation in HPC

In HPC clusters running conda may be a problem.
For those cases there is a a library called 'conda-pack' which allows to run a local lite conda version.
Conda-pack is installed with:

`conda install conda-pack`

When packing a conda environment there must be no package installed in development mode.
If the model workflow is installed in development mode (it can be checked with 'conda list model-workflow') then it must be uninstalled (see previous section). You can install in production mode but it is not required for the conda pack.
Now pack the mwf environment with:

`conda pack -n mwf`

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

`mkdir mwf`

`tar -xzf mwf.tar.gz -C mwf` 

At this point you can access python by runinng:

`mwf/bin/python`

Activate the mwf environment:

`source mwf/bin/activate`

Make a conda unpack

`conda-unpack`

Now copy the whole workflow repository in the remote machine, same repository than before.
Then install it in develop mode with:

`cd workflow`

`../mwf/bin/python setup.py develop`


---

## How to use

print help message:

`mwf --help`

or  `mwf -h`

Run the workflow:

`mwf`

By default, the workflow is run on current directory which is asumed to contain input files: topology, trajectory and the 'inputs.json' file.

A different directory may be specified with the `-dir` or `--working_dir` options.

Input files may be download from an already uploaded project:

`mwf -proj <project id>`

By default, data files will be downloaded from `https://bioexcel-cv19-dev.bsc.es`. Another URL can be specified with the `-url` option.

By default, the whole workflow is run.
The exit may be forced after downloading input files with the `-d` option.
The exit may be forced after downloading input files and/or running mandatory input files processing with the `-s` option.

Run the workflow but include only specific analyses or tools:

`mwf -i <some analysis> <some other analysis> <some tool> ...`

Run the workflow but exlcude some specific analyses:

`mwf -e <some analysis> <some other analysis> ...`

_Other options:_

`-top`: input topology filename. Default: `md.imaged.rot.dry.pdb`

`-traj`: input trajectory filename. Default: `md.imaged.rot.xtc`

`-char`: inputs charges filename. Default: `topology.XXX`

`-in`: inputs filename. Default: `inputs.json`

`-pr`: preprocess protocol. Default: 0

`-trans`: Imaging translation. Default: 0 0 0