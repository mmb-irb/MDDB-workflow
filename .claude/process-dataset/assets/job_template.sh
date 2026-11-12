#!/bin/bash
#SBATCH --partition=mmb_cpu_zen3,mmb_cpu_sphr
#SBATCH --job-name={{DIR}}
#SBATCH --output=mwf_%j.out
#SBATCH --error=mwf_%j.err
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END,FAIL
#SBATCH --time=4:00:00
#SBATCH --mail-user=<your_email>

# `module` is not available on the IRB compute nodes (see references/irb-cluster.md)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate mwf_env

# Run the workflow for this project
# -e: skip the slowest analyses; -lb: local BLAST (see references/launching.md)
mwf run -img -fit -e energies clusters pockets -lb
