## Installation

First clone the workflow repo:

`git clone https://github.com/mmb-irb/MDDB-workflow`

Now create a new environment using the `environment.yml` file in this repo:

`cd MDDB-workflow`<br />
`conda env create --file envs/environment.yml`

Activate the new enviornment

`conda activate mwf_env`

Then install the workflow module in development mode:

`python -m pip install -e ".[dev]"`


## Testing

`pytest --cov-report term --cov=model_workflow test`

## Build

`python -m build --wheel`