#!/bin/bash
# Amber prmtops (complex, receptor, ligand) for one PL-REX complex: plrex_build.sh <complex dir> <work dir>
# Receptor: ff14SB protein + ions + PL-REX extra residues, crystal waters removed (as the energies analysis does).
set -e
cdir=$1; out=$2
FF=$(dirname $(dirname $cdir))/../ff/protein_and_cofactor
FF=$(cd $FF && pwd)
cd $out
[ -f com.prmtop ] && exit 0
# NME caps use old atom names (CH3, HH3x); ff14SB calls them C, H1-H3
grep -E "^(ATOM|HETATM|TER|END)" $cdir/protein.pdb | grep -v " WAT " \
    | sed -E '/ NME /{s/ CH3 NME/ C   NME/; s/HH31 NME/ H1  NME/; s/HH32 NME/ H2  NME/; s/HH33 NME/ H3  NME/}' > receptor_nowat.pdb
cat > tleap.in <<EOT
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
set default PBradii mbondi3
loadamberprep $FF/SEM.prep
NAP = loadmol2 $FF/NAP.mol2
loadamberparams $FF/frcmod.NAP
LIG = loadmol2 lig.mol2
loadamberparams lig.frcmod
rec = loadpdb receptor_nowat.pdb
com = combine {rec LIG}
saveamberparm rec rec.prmtop rec.inpcrd
saveamberparm LIG lig.prmtop lig.inpcrd
saveamberparm com com.prmtop com.inpcrd
quit
EOT
tleap -f tleap.in > tleap.log 2>&1
[ -f com.prmtop ] || { echo "tleap failed: $out"; exit 1; }
