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

1. Merge the master branch into the release branch:

```shell
git checkout release
git merge master
```

2. Wait for the comprehensive test suite to complete. These tests verify workflow integrity and may take a considerable amount of time.

3. Once all tests pass successfully, create and push a new version tag:

```shell
# Create an annotated tag with a version number and message
git tag -a v1.2.3 -m "Release version 1.2.3"

# Push the tag to the remote repository
git push origin v1.2.3
```

4. Create a GitHub release based on this tag through the GitHub web interface, including release notes that detail the changes and improvements.

5. After the release is published, merge the release branch back to master if you made any changes to fix the tests:

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

`make -C docs/ html`

Open `docs/_build/html/index.html` to see the results.

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