# MoDEL Workflow

The aim for this tools is to process raw MD data and obtain standard structure and trajectory files.
These files next to some additional author inputs are analyzed.
Both standard files and analysis results are to be uploaded to the database using the [loader](https://mmb.irbbarcelona.org/gitlab/aluciani/MoDEL-CNS_DB_loader).

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

All you need to start processing your files is a topology file (prmtop, tpr, psf, etc.) and any number of trajectory files (nc, xtc, dcd, etc.). Ideally, every independent trajectory should be in a different directory. Here is an example of a directory we are about to analyze:

raw_topology.prmtop<br />
replica_1/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory.nc<br />
replica_2/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part01.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part02.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part03.nc<br />

Note that here we have 2 independent replicas which share a common topology.<br />
The second replica is actually splitted in 3 consecutive parts but this is not a problem.<br />
The command to start the processing would be as follows:

```bash

mwf -top raw_topology.prmtop -stru raw_topology.prmtop -mdir replica_* -traj *.nc -s

```

The mwf command is the main hanlder of the workflow and its letters stand for 'Model-WorkFlow'.<br />
Now lets explain every argument in this command:<br />
* Argument -top points to the input topology. No big surprises here.<br />
* Argument -stru points to the input structure. Here we do not have any 'pdb' or 'gro' file however. For this reason we also point to the raw topology. We are telling the workflow to get the structure from the topology. To do so it will use some coordinates from the trajectory as well.<br />
* Argument -mdir is telling the workflow which MD directories are to be considered. In this case we want to check all replicas.<br />
* Argument -traj is pointing the workflow towards the input trajectroy files. Note that the path is not relative to our current directory, but to every MD directory.<br />
* Argument -s is just telling the workflow to stop after doing the basic processing steps. Otherwise it would continute running analyses and so on, but we are lacking one last input file to keep going. This is explained further.<br />

If everything is fine after the processing the directory should look like this:

raw_topology.prmtop<br />
topology.prmtop<br />
replica_1/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;structure.pdb<br />
&nbsp;&nbsp;&nbsp;&nbsp;trajectory.xtc<br />
replica_2/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part01.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part02.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part03.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;structure.pdb<br />
&nbsp;&nbsp;&nbsp;&nbsp;trajectory.xtc<br />

Here we must explain a few things more.<br />
First of all, note that trajectory parts have been automatically merged and a structure have been generated for every MD directory.
Note also that the not-raw (or just 'processed') topology is still in amber format (.prmtop) while the processed trajectory is converted to gromacs (.xtc) format. One of the standards in the workflow is that the trajectory is to be called 'trajectory.xtc' as well as the structure is to be called 'structure.pdb'. Thus the trajectory and the structure are to be in 'xtc' and 'pdb' formats respectively. 

<span style="color:blue">The workflow has a wide variety of tools and every tool supports a different range of input formats. Both xtc and pdb formats are extensively used and validated formats. Converting all input formats to a unique format to work along the workflow and thus not having to support every possible structure and trajectory format saves a lot of work.</span>

Conversions between topology formats are difficult however. For this reason the topology is not converted and then used as little as possible. The only processing in topologies is the atom filtering to keep them coherent with both structure and trajectory.

In this example we run the most basic processing, but there are a couple of additional features we may require.
* Filtering: Argument -filt is used to filter atoms away in all files (topology, structure and trajectory). The -filt argument alone applies the default filtering: water and counter ions. However the -filt argument may be followed by some text to apply a custom filtering selection according to [VMD syntax](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html).
* Imaging and fitting: It is not easy to automatize the imaging process. However the workflow is provided with a generic filtering protocol which may be useful in some generic cases. Use the -img argument to image and the -fit argument to fit the trajectory.


---

print help message:

`mwf --help`

or  `mwf -h`

Run the workflow:

`mwf`

By default, the workflow is run on current directory which is asumed to contain input files: topology, trajectory and the 'inputs.json' file.

A different directory may be specified with the `-dir` or `--working_dir` options.

Input files may be download from an already uploaded project:

`mwf -proj <project id>`

By default, data files will be downloaded from `https://mdposit-dev.mddbr.eu`. Another URL can be specified with the `-url` option.

By default, the whole workflow is run.
The exit may be forced after downloading input files with the `-d` option.
The exit may be forced after downloading input files and/or running mandatory input files processing with the `-s` option.

Run the workflow but include only specific analyses or tools:

`mwf -i <some analysis> <some other analysis> <some tool> ...`

Run the workflow but exlcude some specific analyses:

`mwf -e <some analysis> <some other analysis> ...`

_Other options:_

`-stru`: input topology filename. Default: `structure.pdb`

`-traj`: input trajectory filename. Default: `trajectory.xtc`

`-top`: inputs charges filename. Default: `topology.XXX`

`-in`: inputs filename. Default: `inputs.json`

`-img`: is to be imaged. Default: False

`-fit`: is to be fitted. Default: False

`-trans`: Imaging translation. Default: 0 0 0