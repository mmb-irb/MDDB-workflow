# Developing MDDB Workflow

## Installation

First clone the workflow repo and create a new environment using the `environment.yml` file in this repo:

``` shell
git clone https://github.com/mmb-irb/MDDB-workflow mwf
cd mwf
conda env create --file envs/environment.yml
```

Activate the new enviornment and then install the workflow module in development mode:

``` shell
conda activate mwf_env
python -m pip install -e ".[dev,docs]"
```

## Develop on a branch

When working on new features or bug fixes, it's best practice to develop on a separate branch:

```shell
# Create and switch to a new branch
git checkout -b <branch_name>

# Or switch to an existing branch
git switch <branch_name>
```

If you need to apply specific commits from another branch:

```shell
# Apply a specific commit to your current branch
git cherry-pick <commit>
```

Once you're done with your changes, you can create a pull request to merge your branch into the main repository.

## Testing

`pytest --cov-report term --cov=model_workflow test`

## Build wheel

`python -m build --wheel`

## Build docs

`make -C docs/ html`

Open `docs/_build/html/index.html` to see the results.

## Build conda

Test the build with [bioconda-utils](https://bioconda.github.io/contributor/building-locally.html).