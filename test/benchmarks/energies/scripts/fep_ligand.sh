#!/bin/bash
# GAFF2 + AM1-BCC (single-point AM1 on the given pose) for one ligand dir with lig.sdf; writes lig.prmtop/lig.inpcrd.
# Run with mwf_env on PATH: fep_ligand.sh <ligand dir>
set -e
cd $1
[ -f lig.prmtop ] && exit 0
charge=$(python -c "
from rdkit import Chem
m = Chem.MolFromMolFile('lig.sdf', removeHs=False, sanitize=False)
print(sum(a.GetFormalCharge() for a in m.GetAtoms()))")
echo $charge > charge.txt
antechamber -i lig.sdf -fi sdf -o lig.mol2 -fo mol2 -c bcc -nc $charge -at gaff2 -rn LIG -pf y \
    -ek "qm_theory='AM1', maxcyc=0, scfconv=1.d-10, ndiis_attempts=700" > antechamber.log 2>&1
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2
cat > tleap.in <<EOT
source leaprc.gaff2
LIG = loadmol2 lig.mol2
loadamberparams lig.frcmod
saveamberparm LIG lig.prmtop lig.inpcrd
quit
EOT
tleap -f tleap.in > tleap.log 2>&1
[ -f lig.prmtop ] || { echo "tleap failed: $1"; exit 1; }
