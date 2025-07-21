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

The development process is inspired in [trunk-based development](https://trunkbaseddevelopment.com/). When working on new features or bug fixes, it's best practice to develop on a separate branch:

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

Or move a change to another branch:
```shell
# Save your changes to Git's stash
git stash

# Switch to branch
git checkout <branch>

# Apply your stashed changes
git stash apply  
```

Once you're done with your changes, you can create a [pull request (PR)](https://github.com/mmb-irb/MDDB-workflow/pulls) to merge your branch into the main repository. This will trigger the continuous development (CI) tests for some quick check quality checks. It is advised to resolve the PR  with a code review from a peer as a way to check for possible problems, sharing knowledge about the workflow and getting ideas for improvements.

## Making a release

When ready to make a new release, follow these steps:

0. Update the version number in `pyproject.toml`.

1. Merge the master branch into the release branch. Add [skip tests] to the commit message to skip the CI tests:

```shell
git checkout release
git merge master # -m "Merge master into release [skip tests]"
```

2. Wait for the comprehensive test suite to complete. These tests verify workflow integrity and may take a considerable amount of time. If any tests fail, you can fix them in the release branch.

3. Once all tests pass successfully, a new tag is generated automatically based on the `pyproject.toml` version number. This is done by the GitHub Actions workflow defined in [release.yml](.github/workflows/release.yml).

4. After the release is published, merge the release branch back to master if you made any changes to fix the tests:

```shell
git checkout master
git merge release
git push origin master
``` 

## Testing
```shell
# To run a test by its name and optionally a parameter:
pytest test/test_run.py -k test_analysis_execution[pockets]
pytest test/test_run.py -k "TestMWFRun and A01IP and dist"
# To run on a subset of tests use the markers with -m {CI,release}:
pytest -m CI
# To run all tests and generate a coverage report:
pytest --cov-report term --cov=model_workflow -m CI
# Generate the html report and save the console output to report.log while removing color codes
pytest --cov-report xml:coverage/coverage.xml --cov-report html:coverage/ --cov=model_workflow -m release --color=yes | tee >(sed 's/\x1b\[[0-9;]*m//g' > coverage/report.log)
```
```shell
coverage xml -o coverage/coverage.xml
genbadge coverage --name "Coverage" --input-file coverage/coverage.xml  --output-file coverage/coveragebadge.svg
```
## Build wheel

`python -m build --wheel`

## Build docs
Docs are generated automatically for the `master` branch using [Sphinx](https://www.sphinx-doc.org/en/master/). If you have installed the workflow with `python -m pip install -e ".[dev,docs]"`, you can use the following command to build the documentation locally:
```shell
make -C docs/ html
```

Then, you can open `docs/_build/html/index.html` to see the results. This way you can preview the documentation before pushing changes to the repository and without having to wait for the Read the Docs build.

## Build conda

Test the build with [bioconda-utils](https://bioconda.github.io/contributor/building-locally.html):

```shell
conda mambabuild /home/rchaves/repo/biobb/bioconda-recipes/recipes/mddb_workflow -c bioconda
```

Linting:
```shell
cd ~/repo/biobb/bioconda-recipes
bioconda-utils lint  --packages mddb_workflow
```