# MoDEL Workflow

The aim for this tools is to process raw MD data and obtain standard structure and trajectory files.
These files next to some additional author inputs are analyzed.
Both standard files and analysis results are to be uploaded to the database using the [loader](https://mmb.irbbarcelona.org/gitlab/aluciani/MoDEL-CNS_DB_loader).

## Installation

This will install the package and create the entry point `mwf` to be used when executing the workflow.<br />
Note that internet access is required to install the workflow using this protocol. If you are working in a machine where there is not internet then head to the 'Installation in HPC' section.

First clone the workflow repo:

`git clone https://mmb.irbbarcelona.org/gitlab/d.beltran.anadon/MoDEL-workflow.git`

Now create a new environment using the `environment.yml` file in this repo:

`cd MoDEL-workflow`<br />
`conda env create --file environment.yml`

Then install the workflow module in development mode:

`python setup.py develop`


---

## Installation in HPC

In HPC clusters running conda and having internet may be a problem.
For those cases there is a a library called 'conda-pack' which allows to run a local lite conda version.

First install conda-pack in any environment in your locacl machine:

`conda install conda-pack`

Now install in your local machine the conda enviornment as it is explained in the previous section.<br />
Note that there is no need to install the workflow module in this case.

Now pack the mwf environment with:

`conda pack -n mwf --ignore-editable-packages`

This will generate a file called 'mwf.tar.gz'. Copy this file in the remote machine where the workflow must me installed. Then go to that directory and run:

`mkdir mwf_environment`<br />
`tar -xzf mwf.tar.gz -C mwf_environment` 

Activate the mwf environment:

`source mwf_environment/bin/activate`

Make a conda unpack

`conda-unpack`

Now copy the whole workflow repository from your local machine to the remote machine.

`git clone https://mmb.irbbarcelona.org/gitlab/d.beltran.anadon/MoDEL-workflow.git`<br />
`rsync -avP MoDEL-workflow <remote>:<path>`

Then install it in develop mode with:

`cd MoDEL-workflow`<br />
`python setup.py develop`

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

### Processing input files

The command to start the processing would be as follows:

`mwf run -top raw_topology.prmtop -stru raw_topology.prmtop -mdir replica_* -traj *.nc -s`

The 'mwf' command is the main handler of the workflow and its letters stand for 'Model-WorkFlow'.<br />
The 'run' subcommand is the one which triggers the main logic and thus run the workflow.<br />
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
* Imaging and fitting: It is not easy to automatize the imaging process so it is recommended that you manually image your trajectories. However, the workflow is provided with a generic filtering protocol which may be useful in some generic cases. Use the -img argument to image and the -fit argument to fit the trajectory.

Once this process is over some tests are run.<br />
If they all pass then we can continue with the analyses.

### Running the analyses

Before we start, we need one last input file called **inputs.json**.<br />
This file contains burocratic data, MD parameters and some additional metadata which is used by the workflow to adapt the analyses.

In order to generate this file, a Jupyter Notebook to build the file explaining every field in detail is provided. You can find it in the workflow repository, at model_workflow/utils/input_setter.ipynb or open it by simply running the following command:

`mwf inputs`

Fill every field and then run the whole notebook to generate the inputs.json file.
Now we are ready to run the analyses by just running the following command:

`mwf run`

Note that no input files pointers are provided now.<br />
The workflow will guess which are the processed files since they have standard names.<br />
The workflow will also remember which are the MD directories.<br />
In addition, the workflow has some cache in these directories thus remembering some precalculated values and not having to repeat all the process.

Note that this process will also generate some additional files such as 'metadata.json' and 'topology.json' by default. These files are also to be uploaded to the database.

If you want to run only a few specific analyses or exclude some analyses you can use the include (-i) and exclude (-e) arguments.<br />
To see additional arguments and how to use them you can request help by just running `mwf run -h`

Once you are done with this process is time to load your files to the database.<br />
To do so, you must head to the [loader](https://mmb.irbbarcelona.org/gitlab/aluciani/MoDEL-CNS_DB_loader).