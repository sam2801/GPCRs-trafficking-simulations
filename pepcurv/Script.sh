# make peptides
martinize2 -f hecate23.pdb -v -o protein.top -dssp dssp -ff martini3001 -x cg.pdb -elastic
mol=1
rm prot.gro
touch prot.gro
echo 0 | gmx trjconv -f cg.pdb -s cg.pdb -o cg.gro -trans 30 30 10 -quiet

# make bilayer with insane
python2.7 insane.py -f cg.gro -o bilayer.gro -p topol.top -x 60 -y 60 -z 20 -l POPC:60 -l POPS:20 -l CHOL:20 -sol W -salt 0.1 -dm 3 -center

# change topology file
sed -i '9d' topol.top
sed -i -e 's/#include "martini.itp"//g' topol.top
cat << EOF > topol.add
#include "/home/samir/Martini3/martini_v3.0.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_ions_v1.itp"
#include "/home/samir/Martini3/martini_v3.0.0_phospholipids_v1.itp"
#include "/home/samir/Martini3/martini_v3.0_sterols_v1.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_solvents_v1.itp"
#include "molecule_0.itp"
EOF
#mol=$(sed -n '9p' protein.top)
sed -i "7a\ molecule_0    $mol " topol.top
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
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -gpu_id $GPU -nsteps 100
for i in {1..50}
do
echo progress...
echo iteration $i of 50
echo continuing...
gmx grompp -f min.mdp -c min.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -gpu_id $GPU -nsteps 100
done
rm -f \#*
rm -f step*.pdb
for i in {1..20}
do
echo progres...
echo iteration $i of 20
echo continuing...
gmx grompp -f min.mdp -c min.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -gpu_id $GPU -nsteps 300
done
rm -f \#*
rm -f step*.pdb
for i in {1..10}
do
echo progres...
echo iteration $i of 10
echo continuing...
gmx grompp -f min.mdp -c min.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -gpu_id $GPU -nsteps 1000
done
rm -f \#*
rm -f step*.pdb
gmx grompp -f min.mdp -c min.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -gpu_id $GPU

## check energies
echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
xmgrace potential.xvg
rm potential.xvg

# equilibration
cat << EOF > index.input
1 | 2 | 8 | 9 | 10 | 13 | 14 | 15
name 18 LIP
16 | 17
name 19 SOL
q
EOF
gmx make_ndx -f min.gro -quiet < index.input
rm index.input
gmx grompp -f eq.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -n index.ndx -quiet
gmx mdrun -deffnm eq -v -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# simulation
gmx grompp -f md.mdp -p topol.top -c eq.gro -n index.ndx -quiet -o md.tpr
gmx mdrun -v -deffnm md -plumed plumed3.dat -cpi md.cpt -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# vizualize
echo 13 0 | gmx trjconv -f md.xtc -s md.tpr -o md_pbc.xtc -center -pbc whole -dt 200 -quiet
cat << EOF > md_f.pml
load eq.gro
load_traj md_pbc.xtc
remove resname W resname ion
show spheres
show cell
EOF
pymol md_f.pml

# analyze APL of the membrane
fatslim apl -c eq.gro -t md.xtc -n index2.ndx --plot-apl apl.xvg --plot-area area.xvg --hg-group headgroups
fatslim apl -c eq.gro -t md.xtc -n index2.ndx --hg-group headgroups --export-apl-raw md_apl.csv
python3 show_apl_map.py
