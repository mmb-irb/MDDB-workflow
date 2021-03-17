# MoDEL Workflow



## Installation

This will install the package and create the entry point `mwf` to be used when executing the workflow.

It is recommended to create a virtual environment. You can do this using Conda and the `environment.yml` file, with:

```conda env create --file environment.yml```

From the workflow's directory then execute the following:

In DEVELOPMENT:

```python setup.py develop```

To uninstall the develop mode use:

```python setup.py develop -u```

To make notebooks fully operative use:
WARNING: This may break the environment for using the mwf

```conda install -c conda-forge notebook```

In PRODUCTION:

```python setup.py install```

To uninstall the production mode use:

```python setup.py install --record filelist```

```xargs rm -rf < filelist```

---

## Installation in HPC

In HPC clusters running conda may be a problem.
For those cases there is a a library called 'conda-pack' which allows to run a local lite conda version.
Conda-pack is installed with:

```conda install conda-pack```

First of all, you need the environment 'mwf' to be installed in your computer in PRODUCTION mode (see previous section).
Once installed, pack the mwf environment with:

```conda pack -n mwf```

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

```mkdir mwf```

```tar -xzf mwf.tar.gz -C mwf```

At this point you can access python by runinng:

```./mwf/bin/python```

Activate the mwf environment:

```source mwf/bin/activate```

Now copy the whole workflow repository in the remote machine, same repository than before.
Then install it in develop mode with:

```cd workflow```

```../mwf/bin/python setup.py develop```

Finally, some modifications may be required and they must be done by hand

Find the vmd path and modify the vmd executable.

```find ./mwf | grep vmd_LINUX``` (copy the path to the 'lib' directory)

```vim ./mwf/bin/vmd``` (modify the 'defaultvmddir' path as the previous copied path)

At this point the workflow should be operative

---

## How to use

print help message:

```mwf --help```

or  ```mwf -h```

run analysis on current directory, which is asumed to contain input files: topology, trajectory and 'inputs.json'

```mwf```

you can also specify a different directory with the input files using the `-dir` or `--working_dir` options.

download and use input files from already uploaded project

```mwf -p <project id>```

By default the data files will be downloaded from `https://bioexcel-cv19-dev.bsc.es`. Another URL can be specified with the `-url` option.

run single analysis or tool from workflow

```mwf -s <analysis/tool name>```

_Other options:_

`-in`: input filename for the analysis. Default: `inputs.json`

