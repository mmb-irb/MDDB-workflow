# MDDB Workflow

<div align="center">
    <img src="https://github.com/mmb-irb/MDDB-workflow/blob/master/model_workflow/resources/mddb.png?raw=true" width="640" height="320">
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

In HPC clusters running conda and having internet may be a problem.
For those cases there is a a library called 'conda-pack' which allows to run a local lite conda version.

First install conda-pack in any environment in your local machine:

`conda install conda-pack`

Now install in your local machine the conda enviornment as it is explained in the previous section.<br />
Note that there is no need to install the workflow module in this case.

Now pack the mwf environment with:

`conda pack -n mwf --ignore-editable-packages`

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

``` shell
mkdir mwf_env
tar -xzf mwf.tar.gz -C mwf_env
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

## Usage

For detailed instructions on how to use the workflow, please refer to the [usage guide](https://mddb-workflow.readthedocs.io/en/latest/usage.html).

## Development

For instructions about how to develop and test check the [development guide](https://mddb-workflow.readthedocs.io/en/latest/develop.html).
