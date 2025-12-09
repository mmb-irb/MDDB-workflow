# MDDB Workflow

<div align="center">
    <img src="https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/resources/mddb.png?raw=true" width="640" height="320">
</div>

[![Read the Docs (version)](https://img.shields.io/readthedocs/mddb-workflow/latest)](https://mddb-workflow.readthedocs.io/en/latest/)
[![Coverage)](https://raw.githubusercontent.com/mmb-irb/MDDB-workflow/refs/heads/gh-pages/coverage/coveragebadge.svg)](https://mmb-irb.github.io/MDDB-workflow/coverage/)

The aim for this tools is to process raw MD data and obtain standard structure and trajectory files.
These files next to some additional author inputs are analyzed.
Both standard files and analysis results are to be uploaded to the database using the [loader](https://github.com/mmb-irb/MDDB-loader). After that, results can be visualized on the [MDposit](https://mdposit.mddbr.eu/) client.

## Installation

This will install the package and create the entry point `mwf` to be used when executing the workflow.<br />
Note that internet access is required to install the workflow using this protocol. If you are working in a machine where there is not internet then head to the ['Installation in HPC'](#installation-in-hpc) section.

First clone the workflow repo and create a new environment using the `environment.yml` file in this repo:

``` shell
git clone https://github.com/mmb-irb/MDDB-workflow mwf
cd mwf
conda env create --file envs/environment.yml
```

Activate the new enviornment and then install the workflow module in development mode:

``` shell
conda activate mwf_env
python -m pip install -e .
```

At this point your environment is ready to use and it includes the `mwf` command. From now on, you can access the `mwf` command from anywhere in your computer as long as the `mwf_env` environment is activated.


## Installation in HPC

In HPC clusters running conda to install the dependencies cannot be possible in cases where you have no internet access. For this cases you can follow the next steps to install the workflow in your remote machine.

When running the workflow, do not forget to exclude the tasks that need internet accesion by using the flag `-i network`.

### Using conda-pack

Using `conda-pack`, you can create a portable conda environment that can be transferred and deployed on remote machines.

First install conda-pack in any environment in your local machine:

`conda install conda-pack`

Now install in your local machine the conda enviornment as it is explained in the previous section.<br />
Note that there is no need to install the workflow module in this case.

Now pack the mwf environment with:

`conda pack -n mwf_env --ignore-editable-packages`

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

``` shell
rsync -avP mwf_env.tar.gz <remote>:<path>
mkdir mwf_env
tar -xzf mwf_env.tar.gz -C mwf_env
``` 

Activate the mwf environment:

`source mwf_env/bin/activate`

Make a conda unpack

`conda-unpack`

Now copy the whole workflow repository from your local machine to the remote machine.

``` shell
git clone https://github.com/mmb-irb/MDDB-workflow.git mwf
rsync -avP mwf <remote>:<path>
```

Then install it in develop mode with:

``` shell
cd <path>
python -m pip install --no-build-isolation --no-index --no-deps -e .
```

At this point your environment is ready to use and it includes the `mwf` command. From now on, you can access the `mwf` command from anywhere in your computer as long as the `mwf_env` environment is activated.

### Using pixi-pack
If you have `pixi` installed in your local machine, you can use it to create a portable conda environment that can be transferred and deployed on remote machines.
```shell
# Local machine
conda install pixi
pixi init -i envs/environment.yml
pixi lock
pixi-pack --create-executable
rsync -avP environment.sh <remote>:<path>
# Remote machine
./environment.sh
source activate.sh
```

Now copy the whole workflow repository from your local machine to the remote machine.
```shell
# Local machine
git clone https://github.com/mmb-irb/MDDB-workflow.git mwf
rsync -avP mwf <remote>:<path>
# Remote machine
cd <path>
python -m pip install --no-build-isolation --no-index --no-deps -e .
```

### Using containers
Alternatively, you can use containerized versions of the workflow, which bundle all dependencies and the workflow itself into a portable image.

First, pull or build the container image on your local machine, then transfer it to the HPC cluster where you can run the workflow without requiring conda installations.

``` shell
# Using docker
docker save -o mddb_wf.tar ghcr.io/mmb-irb/mddb_wf
scp mddb_wf.tar <remote>:<path>
docker run -v $PWD:/data -w /data ghcr.io/mmb-irb/mddb_wf mwf run

# Using singularity
singularity pull --name mddb_wf.sif docker://ghcr.io/mmb-irb/mddb_wf:latest
scp mddb_wf.sif <remote>:<path>
singularity run -H $PWD -C mddb_wf.sif mwf run
```

## Usage

For detailed instructions on how to use the workflow, please refer to the [usage guide](https://mddb-workflow.readthedocs.io/en/latest/usage.html).

## Development

For instructions about how to develop and test check the [development guide](https://mddb-workflow.readthedocs.io/en/latest/develop.html).

## License

The MDDB Workflow source code is provided under the terms of the Apache
License version 2.0 ([Apache-2.0](LICENSE)).

The MDDB Workflow package (including the workflow with its required
dependencies) is distributed under the terms of the GNU General Public
License version 3 or any later version (GPL-3.0+).
