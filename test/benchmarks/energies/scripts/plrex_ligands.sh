#!/bin/bash
# AM1-BCC charges and GAFF2 parameters for one PL-REX ligand: plrex_ligands.sh <ligand.sdf> <output dir>
set -e
sdf=$1; out=$2
mkdir -p $out && cd $out
[ -f lig.frcmod ] && exit 0
charge=$(awk '/^M  CHG/{for(i=4;i<=NF;i+=2) s+=$(i+1)} END{print s+0}' $sdf)
# Single-point AM1 (maxcyc=0): the PL-REX ligand geometries are already optimized in the complex
antechamber -i $sdf -fi sdf -o lig.mol2 -fo mol2 -c bcc -nc $charge -at gaff2 -rn LIG -pf y \
    -ek "qm_theory='AM1', maxcyc=0, scfconv=1.d-10, ndiis_attempts=700" > antechamber.log 2>&1
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2
echo $charge > charge.txt
