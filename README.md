# MoDEL Workflow



## How to make the workflow executable

change permissions for `model_wf.py` file by doing

```chmod +x model_wf.py```

then, create a symlink

```ln -s model_wf.py model_wf```

finally, add directory that contains the symlink to `PATH`

```export PATH=$(pwd):$PATH```

and you're done!

---

## How to use

print help message:

```model_wf --help```

**actions**

* run analysis on directory which contains topology, trajectory and input files

```model_wf -analyze_directory <directory path>```

you can also specify the trajectory filename with `-topology_filename` and the topology filename with `-topology_filename` (they must be located in the directory path specified).

* run analysis from already uploaded project

```model_wf -analyze_project <project name> <directory name> <URL (optional)>```

you can again specify trajectory and topology filenames.

* run single analysis from workflow

```model_wf -single_analysis <analysis name>```

you can specify the directory with all the needed files with `-working_dir`, otherwise it will use the current directory. Topology and trajectory files can be specified like previously mentioned. 

Analysis name choices: `{rmsd,rmsf,rgyr,pca,pca_contacts,rmsd_per_residue,rmsd_pairwise,distance_per_residue,hydrogen_bonds,energies,pockets}`

_Other options:_

`-output_analysis_filename`: output filename for the analysis.

`-pca_eigenval_file/-pca_eigenvec_file`: eigenvalue and eigenvector filenames for PCA analysis.

`--analysis_prep`: if used, before the specified analysis it will run preparation where it generates metadata, computes interfaces, etc. 

`-interface_cutoff_distance`: cut-off distance to be used when computing interfaces.

`-inputs_filename`: input data filename.

`-metadata_filename`: metadata filename.

