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

First of all, you need the environment 'mwf' to be installed in your computer in PRODUCTION mode (see previous section).
Once installed, pack the mwf environment with:

`conda pack -n mwf`

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

`mkdir mwf`

`tar -xzf mwf.tar.gz -C mwf` 

At this point you can access python by runinng:

`mwf/bin/python`

Activate the mwf environment:

`source mwf/bin/activate`

Now copy the whole workflow repository in the remote machine, same repository than before.
Then install it in develop mode with:

`cd workflow`

`../mwf/bin/python setup.py develop`

### Check everything works fine

Sometimes, some tools in the packed environment may have wrong paths to their executables

- Check VMD

Activate the mwf enviornment and type `vmd`

If an error similar to "/slgpfs/projects/bsc23/bsc23088/mwf/bin/vmd: line 505: /opt/anaconda1anaconda2anaconda3/lib/vmd_LINUX: No existe el fichero o el directorio" is return then do the following:

Find the vmd path and modify the vmd executable.

`find mwf | grep vmd_LINUX` (copy the path to the 'lib' directory, not to the vmd executable)

`vim mwf/bin/vmd` (modify the 'defaultvmddir' path as the previous copied path, relative to mwf)

At this point the workflow should be operative

- Check Gromacs

Activate the mwf enviornment and type `gmx`

If an error similar to "/opt/conda/conda-bld/gromacs_1618496378051/_h_env_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold_place/bin.SSE2/GMXRC: No existe el fichero o el directorio" is return then do the following:

`find "$(pwd -P)" | grep GMXRC` (copy the path which contains all the bin."whatever" directories)

`vim mwf/bin/gmx` (modify all wrong paths as the previous copied path, absolute)

You may also have to modify the atom masses file to add elements typed with both upper case letters (e.g. 'ZN' in addition to 'Zn')

`vim mwf/share/gromacs/top/atommass.dat`

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

`mwf -p <project id>`

By default, data files will be downloaded from `https://bioexcel-cv19-dev.bsc.es`. Another URL can be specified with the `-url` option.

By default, the workflow is run when all required stuff is already downloaded. The exit may be forced before running the workflow with the `-s` option.

Run the workflow but include only specific analyses or tools:

`mwf -i <some analysis> <some other analysis> <some tool> ...`

Run the workflow but exlcude some specific analyses:

`mwf -e <some analysis> <some other analysis> ...`

_Other options:_

`-top`: input topology filename. Default: `md.imaged.rot.dry.pdb`

`-traj`: inputs filename. Default: `md.imaged.rot.xtc`

`-in`: inputs filename. Default: `inputs.json`

`-pr`: preprocess protocol. Default: 0

