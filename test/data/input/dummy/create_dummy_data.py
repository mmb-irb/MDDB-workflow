# Script to create dummy molecular dynamics data files for GROMACS and AMBER simulations.
# Not meant to be run, just here for provenance.
import os
import subprocess
import numpy as np
import MDAnalysis as mda


def create_amber_files(output_dir, water=False):
    """ Creates AMBER files (prmtop, inpcrd, pdb) for an Ala-Ala dipeptide in water. """

    leap_script_path = os.path.join(output_dir, 'leap.in')
    with open(leap_script_path, 'w') as f:
        f.write("source leaprc.protein.ff14SB\n\n")
        f.write("mol = sequence { NALA ALA CALA }\n\n")
        if water:
            f.write("source leaprc.water.tip3p\n")
            f.write("solvatebox mol TIP3PBOX 10.0 iso\n\n")
            f.write("addions mol Na+ 0\n")
            f.write("addions mol Cl- 0\n\n")
        f.write("saveamberparm mol ala_ala.prmtop ala_ala.inpcrd\n\n")
        f.write("quit\n")
    print(f"Created tleap script: {leap_script_path}")

    prmtop_path = os.path.join(output_dir, 'ala_ala.prmtop')
    inpcrd_path = os.path.join(output_dir, 'ala_ala.inpcrd')

    try:
        # Run tleap
        subprocess.run(
            ['tleap', '-f', 'leap.in'],
            check=True, capture_output=True, text=True, cwd=output_dir
        )
        print(f"Created AMBER prmtop file: {prmtop_path}")
        print(f"Created AMBER inpcrd file: {inpcrd_path}")

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        if isinstance(e, subprocess.CalledProcessError):
            print(f"tleap stdout:\n{e.stdout}")
            print(f"tleap stderr:\n{e.stderr}")
        return None

    # Extract PDB file using ambpdb
    pdb_path = 'ala_ala.pdb'

    with open(inpcrd_path, 'r') as inpcrd_file:
        result = subprocess.run(
            ['ambpdb', '-p', 'ala_ala.prmtop'],
            stdin=inpcrd_file,
            capture_output=True,
            text=True,
            check=True,
            cwd=output_dir
        )

    with open(pdb_path, 'w') as pdb_file:
        pdb_file.write(result.stdout)

    print(f"Created PDB file: {pdb_path}")


    # Create a simple NetCDF trajectory with 3 frames
    nc_path = os.path.join(output_dir, 'trajectory.nc')
    u = mda.Universe(prmtop_path, inpcrd_path)
    with mda.Writer(nc_path, n_atoms=u.atoms.n_atoms) as W:
        # First frame (t=0) - original positions
        W.write(u.atoms)

        # Second frame (t=1) - slightly perturbed positions
        u.atoms.positions += np.random.randn(u.atoms.n_atoms, 3) * 0.1
        W.write(u.atoms)

        # Third frame (t=2) - more perturbation
        u.atoms.positions += np.random.randn(u.atoms.n_atoms, 3) * 0.1
        W.write(u.atoms)
        print(f"Created AMBER NetCDF trajectory: {nc_path}")


    # Clean up intermediate files
    os.remove(leap_script_path)
    for temp_file in ['leap.log']:
        temp_path = os.path.join(output_dir, temp_file)
        if os.path.exists(temp_path):
            os.remove(temp_path)

    return pdb_path


def create_gromacs_files(pdb_path, output_dir):
    """ Creates GROMACS files (gro, top) from AMBER PDB using gmx pdb2gmx. """

    gro_path = os.path.join(output_dir, 'ala_ala.gro')
    top_path = os.path.join(output_dir, 'topol.top')

    try:
        # Run gmx pdb2gmx
        # Use force field 6 (AMBER99SB-ILDN) and water model 1 (TIP3P) as defaults
        # The -ignh flag ignores hydrogen atoms in the input (will be added by pdb2gmx)
        subprocess.run(
            ['gmx', 'pdb2gmx', '-f', pdb_path, '-o', gro_path, '-p', top_path,
             '-water', 'tip3p', '-ff', 'amber99sb-ildn', '-ignh'],
            check=True, capture_output=True, text=True,
            input='\n'  # Provide empty input in case it asks for selections
        )
        print(f"Created GROMACS gro file: {gro_path}")
        print(f"Created GROMACS top file: {top_path}")

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        if isinstance(e, subprocess.CalledProcessError):
            print(f"gmx pdb2gmx stderr:\n{e.stderr}")
        return

    # Create a simple XTC trajectory with 3 frames
    xtc_path = os.path.join(output_dir, 'trajectory.xtc')
    try:
        u = mda.Universe(gro_path)
        with mda.Writer(xtc_path, n_atoms=u.atoms.n_atoms) as W:
            # First frame (t=0) - original positions
            W.write(u.atoms)

            # Second frame (t=1) - slightly perturbed positions
            u.atoms.positions += np.random.randn(u.atoms.n_atoms, 3) * 0.1
            W.write(u.atoms)

            # Third frame (t=2) - more perturbation
            u.atoms.positions += np.random.randn(u.atoms.n_atoms, 3) * 0.1
            W.write(u.atoms)
        print(f"Created GROMACS XTC trajectory: {xtc_path}")
    except Exception as e:
        print(f"Error creating XTC trajectory: {e}")

    # Create a simple TPR file (optional, requires mdp file)
    mdp_path = os.path.join(output_dir, 'grompp.mdp')
    mdp_out = os.path.join(output_dir, 'mdout.mdp')
    tpr_path = os.path.join(output_dir, 'ala_ala.tpr')

    with open(mdp_path, 'w') as f:
        f.write("integrator  = md\n")
        f.write("nsteps      = 10\n")
        f.write("dt          = 0.001\n")
        f.write("cutoff-scheme = Verlet\n")
        f.write("rlist       = 0.1\n")
        f.write("rvdw        = 0.1\n")
        f.write("rcoulomb    = 0.1\n")
        f.write("pbc         = xyz\n")

    try:
        os.remove
        subprocess.run(
            ['gmx', 'grompp', '-f', 'grompp.mdp', '-c', 'ala_ala.gro', '-p', 'topol.top',
             '-o', 'ala_ala.tpr', '-maxwarn', '10', '-po', 'mdout.mdp'],
            check=True, capture_output=True, text=True, cwd=output_dir
        )
        print(f"Created GROMACS TPR file: {tpr_path}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Warning: Could not create TPR file.")
        if isinstance(e, subprocess.CalledProcessError):
            print(f"gmx grompp stderr:\n{e.stderr}")

    # Clean up temporary files
    for file in os.listdir(output_dir):
        if file.endswith('mdp') or file.startswith('#'):
            os.remove(os.path.join(output_dir, file))
    os.remove('posre.itp')



if __name__ == "__main__":
    # Create output directories if they don't exist
    gromacs_output_dir = "gromacs"
    amber_output_dir = "amber"
    os.makedirs(gromacs_output_dir, exist_ok=True)
    os.makedirs(amber_output_dir, exist_ok=True)

    print("--- Creating AMBER dummy files ---")
    pdb_path = create_amber_files(amber_output_dir)

    print("\n--- Creating GROMACS dummy files ---")
    create_gromacs_files(pdb_path, gromacs_output_dir)

    print("\nDummy data creation complete.")
