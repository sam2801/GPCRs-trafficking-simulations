# make bilayer with insane
python2.7 insane.py -o bilayer.gro -p topol.top -x 60 -y 60 -z 20 -hole 3 -l POPC:60 -l POPS:20 -l CHOL:20 -sol W -salt 0.1

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
2 | 3 | 4
name 7 LIP
5 | 6
name 8 SOL
q
EOF
gmx make_ndx -f min.gro -quiet < index.input
rm index.input
gmx grompp -f eq.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -n index.ndx -quiet -maxwarn 1
gmx mdrun -deffnm eq -v -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# encurv
cat << EOF > plumed3.dat
# Define dummy atom at the center of curvature
piv: FIXEDATOM AT=45,45,12

# Atoms within 13nm of the dummy atom in the membrane
membr: GROUP NDX_FILE=index2.ndx NDX_GROUP=patch

# Restraint groups
#r1: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_1
#r2: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_2
#r3: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_3
#r4: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_4
#r5: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_5
#r6: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_6
#r7: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_7
#r8: GROUP NDX_FILE=index2.ndx NDX_GROUP=restrain_8

# Center of mass for restraint squares
#c1: COM ATOMS=r1
#c2: COM ATOMS=r2
#c3: COM ATOMS=r3
#c4: COM ATOMS=r4
#c5: COM ATOMS=r5
#c6: COM ATOMS=r6
#c7: COM ATOMS=r7
#c8: COM ATOMS=r8

# Distance between centers of mass and virtual atom
#d1: DISTANCE ATOMS=piv,c1
#d2: DISTANCE ATOMS=piv,c2
#d3: DISTANCE ATOMS=piv,c3
#d4: DISTANCE ATOMS=piv,c4
#d5: DISTANCE ATOMS=piv,c5
#d6: DISTANCE ATOMS=piv,c6
#d7: DISTANCE ATOMS=piv,c7
#d8: DISTANCE ATOMS=piv,c8

# Define upper wall
#uwal1: UPPER_WALLS ARG=d1 AT=6 KAPPA=150
#uwal2: UPPER_WALLS ARG=d2 AT=6 KAPPA=150
#uwal3: UPPER_WALLS ARG=d3 AT=6 KAPPA=150
#uwal4: UPPER_WALLS ARG=d4 AT=6 KAPPA=150
#uwal5: UPPER_WALLS ARG=d5 AT=6 KAPPA=150
#uwal6: UPPER_WALLS ARG=d6 AT=6 KAPPA=150
#uwal7: UPPER_WALLS ARG=d7 AT=6 KAPPA=150
#uwal8: UPPER_WALLS ARG=d8 AT=6 KAPPA=150

# Define lower wall
#lwal1: LOWER_WALLS ARG=d1 AT=4 KAPPA=150
#lwal2: LOWER_WALLS ARG=d2 AT=4 KAPPA=150
#lwal3: LOWER_WALLS ARG=d3 AT=4 KAPPA=150
#lwal4: LOWER_WALLS ARG=d4 AT=4 KAPPA=150
#lwal5: LOWER_WALLS ARG=d5 AT=4 KAPPA=150
#lwal6: LOWER_WALLS ARG=d6 AT=4 KAPPA=150
#lwal7: LOWER_WALLS ARG=d7 AT=4 KAPPA=150
#lwal8: LOWER_WALLS ARG=d8 AT=4 KAPPA=150

# Define EnCurv collective variable with R=8
sec1: ENCURV ATOMS=piv,membr NBINS=50 R=8

# Harmonic restraint on radial component
restr1: RESTRAINT ARG=sec1.val KAPPA=100 AT=8 STRIDE=2

# Harmonic restraint on angular component
#angrestr: RESTRAINT ARG=sec1.angle KAPPA=10000 AT=0 STRIDE=2

# Print RMSD
PRINT STRIDE=100 ARG=sec1.rmsd,sec1.angle FILE=COLVAR
EOF
# echo 0 | gmx trjconv -f eq.xtc -s eq.tpr -fr frame.ndx -o frame.gro -quiet
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on select.ndx -select 'same residue as (resname POPC or resname POPS or resname CHOL) and (within 10 of [45, 45, 10] or within 10 of [45, 45, 7.5] or within 10 of [45, 45, 12.5])'
sed -i -e 's/.*same_residue_as_(resname_POPC_or_resname_POPS_or_resname_CHOL)_and.*/[ patch ]/g' select.ndx
count=1
for i in 0 27 54
do
for j in 0 27 54
do
a=$((i+6))
b=$((j+6))
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on select$count.ndx -select "same residue as (resname POPC or resname POPS or resname CHOL) and x > ${i} and x < ${a} and y > ${j} and y < ${b}"
count=$((count+1))
done
done
rm select5.ndx
square=1
for i in 1 2 3 4 6 7 8 9
do
sed -i -e "s/.*same_residue_as_(resname_POPC_or_resname_POPS_or_resname_CHOL)_and.*/[ restrain_${square} ]/g" select$i.ndx
square=$((square+1))
done
cat index.ndx select.ndx select1.ndx select2.ndx select3.ndx select4.ndx select6.ndx select7.ndx select8.ndx select9.ndx > index2.ndx
rm \select*
gmx grompp -f md.mdp -p topol.top -c eq.gro -r eq.gro -n index2.ndx -quiet -o md.tpr -maxwarn 1
i=1
until gmx mdrun -v -deffnm md -plumed plumed3.dat -cpi md.cpt -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet; do
    echo "Simulation crashed ${i} times."
    rm -f step*.pdb
    rm -f \#*
    i=$((i+1))
    sleep 1
done

tail -n +2 COLVAR > colvar_graph.txt
python3 colvar_graph.py

# comparison
cat << EOF > plumed_2d.dat
# Define dummy atom at the center of curvature
piv: FIXEDATOM AT=30,30,13

# Atoms within 13nm of the dummy atom in the membrane
membr: GROUP NDX_FILE=index2.ndx NDX_GROUP=patch

# Define EnCurv collective variable with R=8
sec1: ENCURV ATOMS=piv,membr NBINS=50 R=8

# Harmonic restraint on radial component
restr1: RESTRAINT ARG=sec1.val KAPPA=100 AT=8 STRIDE=2

# Harmonic restraint on angular component
#angrestr: RESTRAINT ARG=sec1.angle KAPPA=10000 AT=0 STRIDE=2

# Print RMSD
PRINT STRIDE=100 ARG=sec1.rmsd,sec1.angle FILE=COLVAR
EOF
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on select.ndx -select 'same residue as (resname POPC or resname POPS or resname CHOL) and (within 8 of [30, 30, 10] or within 8 of [30, 30, 7.5] or within 8 of [30, 30, 12.5])'

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
fatslim apl -c eq.gro -t md.xtc -n index2.ndx --plot-apl apl.xvg --plot-area area.xvg --hg-group headgroups
