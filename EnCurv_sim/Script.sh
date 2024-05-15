# make bilayer with insane
python2.7 insane.py -o bilayer.gro -p topol.top -x 60 -y 60 -z 20 -l POPC:60 -l POPS:20 -l CHOL:20 -sol W -salt 0.1

# change topology file
sed -i -e 's/#include "martini.itp"//g' topol.top
cat << EOF > topol.add
#include "/home/samir/Martini3/martini_v3.0.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_ions_v1.itp"
#include "/home/samir/Martini3/martini_v3.0.0_phospholipids_v1.itp"
#include "/home/samir/Martini3/martini_v3.0_sterols_v1.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_solvents_v1.itp"
EOF
cat topol.add topol.top > tmp
rm topol.add
mv tmp topol.top
sed -i -e 's/NA+/NA /g' topol.top
sed -i -e 's/CL-/CL /g' topol.top

# change structure file
sed -i -e 's/NA+ /NA  /g' bilayer.gro
sed -i -e 's/ NA+/  NA/g' bilayer.gro
sed -i -e 's/CL- /CL  /g' bilayer.gro
sed -i -e 's/ CL-/  CL/g' bilayer.gro

# minimization
gmx grompp -f min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet

# equilibration
cat << EOF > index.input
2 | 3 | 4 | 9 | 10 | 11
name 15 LIP
5 | 6 | 7 | 8 | 12 | 13 | 14
name 16 SOL
q
EOF
gmx make_ndx -f bilayer.gro -quiet < index.input
rm index.input
gmx grompp -f eq.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -n index.ndx -quiet
gmx mdrun -deffnm eq -v -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# encurv
cat << EOF > plumed3.dat
# Define dummy atom at the center of curvature
piv: FIXEDATOM AT=30,30,13

# Atoms within 13nm of the dummy atom in the membrane
membr: GROUP NDX_FILE=index2.ndx NDX_GROUP=patch

# Define EnCurv collective variable with R=10
sec1: ENCURV ATOMS=piv,membr NBINS=50 R=10

# Harmonic restraint on radial component
restr1: RESTRAINT ARG=sec1.val KAPPA=100 AT=10 STRIDE=2

# Harmonic restraint on angular component
#angrestr: RESTRAINT ARG=sec1.angle KAPPA=10000 AT=0 STRIDE=2

# Print RMSD
PRINT STRIDE=100 ARG=sec1.rmsd,sec1.angle FILE=COLVAR
EOF
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on select.ndx -select 'same residue as (resname POPC or resname POPS or resname CHOL) and (within 8 of [30, 30, 10] or within 8 of [30, 30, 7.5] or within 8 of [30, 30, 12.5])'
sed -i -e 's/.*same_residue_as.*/[ patch ]/g' select.ndx
cat index.ndx select.ndx > index2.ndx
rm \select*
gmx grompp -f md.mdp -p topol.top -c eq.gro -n index2.ndx -quiet -o md.tpr
gmx mdrun -v -deffnm md -plumed plumed3.dat -cpi md.cpt -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

sed -e 's/.*FIELDS.*/ /g' COLVAR > colvar_graph.txt
python3 colvar_graph.py


# vizualize
echo 2 0 | gmx trjconv -f md.xtc -s md.tpr -o md_pbc.xtc -center -pbc whole -dt 200 -quiet
cat << EOF > md.pml
load md.gro
load_traj md_pbc.xtc
remove resname W resname WF
show spheres
EOF
pymol md.pml

# analyze APL of the membrane
fatslim apl -c eq.gro -t md.xtc -n index3.ndx --plot-apl apl.xvg --plot-area area.xvg --hg-group headgroups
fatslim apl -c eq.gro -t md.xtc -n index3.ndx --hg-group headgroups --export-apl-raw md_apl.csv
python3 show_apl_map.py
