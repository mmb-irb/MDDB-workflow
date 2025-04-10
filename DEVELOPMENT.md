## Installation

First clone the workflow repo:

`git clone https://github.com/mmb-irb/MDDB-workflow mwf`

Now create a new environment using the `environment.yml` file in this repo:

``` shell
cd MDDB-workflow
conda env create --file envs/environment.yml
```

Activate the new enviornment and then install the workflow module in development mode:

``` shell
conda activate mwf_env
python -m pip install -e ".[dev]"
```

## Testing

`pytest --cov-report term --cov=mddb_wf test`

## Build

`python -m build --wheel`