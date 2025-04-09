# MDDB Workflow

The aim for this tools is to process raw MD data and obtain standard structure and trajectory files.
These files next to some additional author inputs are analyzed.
Both standard files and analysis results are to be uploaded to the database using the [loader](https://github.com/mmb-irb/MDDB-loader).

## Installation

This will install the package and create the entry point `mwf` to be used when executing the workflow.<br />
Note that internet access is required to install the workflow using this protocol. If you are working in a machine where there is not internet then head to the 'Installation in HPC' section.

`conda install bioconda::mddb`

At this point your environment is ready to use and it includes the `mwf` command.<br />
From now on, you can access the `mwf` command from anywhere in your computer as long as the `mwf_env` environment is activated.

---

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

`mkdir mwf_env`<br />
`tar -xzf mwf.tar.gz -C mwf_env` 

Activate the mwf environment:

`source mwf_env/bin/activate`

Make a conda unpack

`conda-unpack`

Now copy the whole workflow repository from your local machine to the remote machine.

`git clone https://mmb.irbbarcelona.org/gitlab/d.beltran.anadon/MoDEL-workflow.git`<br />
`rsync -avP MoDEL-workflow <remote>:<path>`

Then install it in develop mode with:

`cd MoDEL-workflow`<br />
`python setup.py develop`

At this point your environment is ready to use and it includes the `mwf` command.<br />
From now on, you can access the `mwf` command from anywhere in your computer as long as the `mwf_env` environment is activated.

---

## How to use

All you need to start processing your files is a topology file (prmtop, tpr, psf, etc.) and any number of trajectory files (nc, xtc, dcd, etc.). Every independent trajectory must be in a different directory, where several output files will be generated for every trajectory. Here is an example of a directory we are about to analyze:

raw_topology.prmtop<br />
replica_1/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory.nc<br />
replica_2/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part01.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part02.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part03.nc<br />

Note that here we have 2 independent replicas which share a common topology.<br />
The second replica is actually splitted in 3 consecutive parts but this is not a problem.<br />

> :warning: WARNING<br />
> Keeping data organized like in the example is important for the workflow to work. In order to keep your trajectories safe we recommend making a copy or using symlinks as workflow inputs. e.g. `ln -s ~/my_data/my_trajectory.nc replica_1/raw_trajectory.nc`.
Note that the workflow will generate output files where your input trajectories are.

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

If everything is fine after the processing some <span style="color: green;">new files</span> may have been generated and now the directory should look like this:

raw_topology.prmtop<br />
<span style="color: green;">topology.prmtop</span><br />
replica_1/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;<span style="color: green;">structure.pdb</span><br />
&nbsp;&nbsp;&nbsp;&nbsp;<span style="color: green;">trajectory.xtc</span><br />
replica_2/<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part01.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part02.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;raw_trajectory_part03.nc<br />
&nbsp;&nbsp;&nbsp;&nbsp;<span style="color: green;">structure.pdb</span><br />
&nbsp;&nbsp;&nbsp;&nbsp;<span style="color: green;">trajectory.xtc</span><br />

Here we must explain a few things more.<br />
First of all, note that trajectory parts have been automatically merged and a structure have been generated for every MD directory.
Note also that the not-raw (or just 'processed') topology is still in amber format (.prmtop) while the processed trajectory is converted to gromacs (.xtc) format. One of the standards in the workflow is that the trajectory is to be called 'trajectory.xtc' as well as the structure is to be called 'structure.pdb'. Thus the trajectory and the structure are to be in 'xtc' and 'pdb' formats respectively.

> The workflow has a wide variety of tools and every tool supports a different range of input formats. Both xtc and pdb formats are extensively used and validated formats. Converting all input formats to a unique format to work along the workflow and thus not having to support every possible structure and trajectory format saves a lot of work.<br />
>
> Conversions between topology formats are difficult however. For this reason the topology is not converted and then used as little as possible. The only processing in topologies is the atom filtering to keep them coherent with both structure and trajectory.

In this example we run the most basic processing, but there are a couple of additional features we may require.
* Filtering: Argument -filt is used to filter atoms you want to keep in all files (topology, structure and trajectory). The -filt argument alone applies the default filtering: water and counter ions. However the -filt argument may be followed by some text to apply a custom filtering selection according to [VMD syntax](https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html). Here is an example on how to filter away everything which is not protein or nucleic acids:<br />
`mwf run (...) -filt 'protein or nucleic'`
* Imaging and fitting: It is not easy to automatize the imaging process so it is recommended that you manually image your trajectories. However, the workflow is provided with a generic imaging protocol which may be useful in some generic cases. Use the -img argument to image and the -fit argument to fit the trajectory.

Once this process is over some tests are run.<br />
If they all pass then we can continue with the analyses.

### Running the analyses

Before we start, we need the **inputs file**.<br />
This file contains burocratic data, MD parameters and some additional metadata which is used by the workflow to adapt the analyses.

In order to generate this file, a template to build the file explaining every field in detail is provided. You can find it in the workflow repository, at [model_workflow/resources/inputs_file_template.yml](model_workflow/resources/inputs_file_template.yml) or open it by simply running the following command:

`mwf inputs`

Fill every field and then run the whole notebook to generate the inputs file.
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

### Tests and other checking processes

Some tests and checkings are run along the workflow to ensure data quality.<br />
By default, a falied test will kill the workflow in the spot.<br />
Tests may be skipped with the -t or --trust argument.<br />
Tests may be allowed to fail with the -m or --mercy argument.<br />
Either if a test is skipped or failed it will write a warning log in the metadata.<br />
Note that some tests may not be skipped, but only allowed to fail and vice versa.<br />
Both the trust and the mercy arguments have the same behaviour:<br />
They may be followed by the name of the test to be skipped/allowed to fail.<br />
If they are passed alone then they are applied to all available tests.<br />
For example, if you want to allow bonds to be wrong, you may skip both the stable bonds test and the coherent bonds test:
`mwf run --trust stabonds cohbonds`

These are the available tests:<br />

- **stabonds** - Stable bonds test

Check atom bonds to be stable along the trajectory.<br />
Bonds are mined from the topology when possible.<br />
Otherwise they are guessed by atom distance and radii ([VMD logic](https://www.ks.uiuc.edu/Research/vmd/current/ug/node27.html). under the hood) 10 frames along the trajectory are checked and bonds happening in most of these frames are the ones to be considered real. This way we avoid taking atom clashes as real bonds.<br />
Then, either if bonds are mined or guessed, we search for the first frame in the trajectory which respects all these bonds using again the atom distance and radii logic. This frame is set as the "reference frame" and it is usually found among the first 10 frames in the trajectory. If no frame matches all bonds after checking the first 100 frames then the test fails.<br />
If this test is skipped or allowed to fail then the first frame in the trajectory is set as the reference frame.

- **cohbonds** - Coherent bonds

Check number of bonds per atom to be coherent according to chemistry knowledge:<br />
<li>Hydrogen atoms must have 1 and only 1 bond</li>
<li>Oxygen atoms must have between 1 and 2 bonds</li>
<li>Nitrogen atoms must have between 1 and 4 bonds</li>
<li>Carbon atoms must have between 2 and 4 bonds</li>
<li>Sulfur atoms must have between 1 and 6 bonds</li>
<li>Phosphorus atoms must have between 2 and 6 bonds</li>

If any of these rules is not respected then the test fails.

- **intrajrity** - Trajectory integrity

Make sure there are no sudden jumps in the middle of the trajectory due to imaging problems.<br />
Compute the RMSD between every pair of consecuitve frames in the trajectory.<br />
Then calculate the standard deviation among all RMSD differences.<br />
If there is at least one jump which is greater than 9 times the standard deviation then the test fails.<br />
First frames are excepcionally allowed to reach this limit since they may be part of the equilibration proccess.<br />

- **elements** - Correct elements (skip only)

Set wrong or missing atom elements to make them standard.<br />
This process relies in guessing mostly but it is smart enough to not set alpha carbons as calcium (Ca).<br />
This process logs warnings for every element which is guessed to be wrong.<br />
If this process is skipped it raises warning anyway but it keeps the original elements intact.<br />

- **refseq** - Reference sequence match (mercy only)

Make sure all protein chains are matched with their corresponding UniProt id.<br />
UniProt ids may be passed through the inputs file.<br />
Additional UniProt ids may be mined from the PDB in case PDB ids are passed through the inputs file.<br />
If there is a protein chain which finds no match among the available UniProt ids then a BLAST is run against reviewed UniProt sequences only.<br />
If the BLAST also fails then the matching process fails.<br />
If this process is allowed to fail then the unmatched protein region will remain with no UniProt reference.<br />
Note that proteins which are not to match anything such as antibodies or artifical constructs are to be tagged as 'no referable' in the inputs file.

- **interact** - Stable interactions (mercy only)

Make sure defined interactions are actually happening and stable enough to be computed.<br />
In order to find interface residues, interactions are checked to happend 100 frames along the trajectory.<br />
Two residues are considered to be in the interface when they are close enough at least in one of these frames.<br />
If no interface residues are found then the process is killed and this can not be allowed.<br />
This means you defined an interaction which does not exist and thus it must be removed from the inputs file.<br />
However, it may happen that an interaction has interface residues but it is almost not happening.<br />
For example, a ligand which is fitted in a protein pocket but it leaves its place as soon as the trajectroy starts playing.<br />
For this reason, it is computed the number of frames that the interaction actually happens.<br />
If the interaction takes place less than the 10% of the total trajectory then the process fails.<br />
If it is allowed to fail then the wrong interaction is removed from metadata anyway and interaction analyses are not run for this specific interaction.