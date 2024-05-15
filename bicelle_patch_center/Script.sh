# make bilayer with insane
# try DCPC
python2.7 insane.py -o bilayer.gro -p topol.top -x 70 -y 70 -z 20 -l POPC:60 -l POPS:20 -l CHOL:20 -sol W -salt 0.1

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
gmx mdrun -deffnm min -v -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet

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

# isolate circular patch
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on pc.ndx -select 'same residue as resname POPC and (within 30 of [35, 35, 10] or within 30 of [35, 35, 7.5] or within 30 of [35, 35, 12.5])'
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on ps.ndx -select 'same residue as resname POPS and (within 30 of [35, 35, 10] or within 30 of [35, 35, 7.5] or within 30 of [35, 35, 12.5])'
gmx select -f eq.gro -selrpos res_com -s eq.tpr -quiet -on col.ndx -select 'same residue as resname CHOL and (within 30 of [35, 35, 10] or within 30 of [35, 35, 7.5] or within 30 of [35, 35, 12.5])'
sed -i -e 's/.*same_residue_as_resname.*/[ POPC_bicelle ]/g' pc.ndx
sed -i -e 's/.*same_residue_as_resname.*/[ POPS_bicelle ]/g' ps.ndx
sed -i -e 's/.*same_residue_as_resname.*/[ CHOL_bicelle ]/g' col.ndx
cat << EOF > sys.ndx
[ Bicelle ]
EOF
tail -n +2 pc.ndx >> sys.ndx
tail -n +2 ps.ndx >> sys.ndx
tail -n +2 col.ndx >> sys.ndx
mol_POPC=$(($(ls | tail -n +2 pc.ndx | wc -w)/12))
mol_POPS=$(($(ls | tail -n +2 ps.ndx | wc -w)/12))
mol_CHOL=$(($(ls | tail -n +2 col.ndx | wc -w)/9))
cat sys.ndx pc.ndx ps.ndx col.ndx > bc.ndx
rm pc.ndx
rm ps.ndx
rm col.ndx
rm sys.ndx
echo 0 | gmx trjconv -f eq.gro -s eq.tpr -n bc.ndx -o bc.gro -quiet

# solvate
python2 insane.py -f bc.gro -o sol.gro -p sol.top -x 68.14030 -y 68.14030 -z 19.91084 -sol W -salt 0.1

# topology
sed -i -e 's/#include "martini.itp"//g' sol.top
cat << EOF > topol.add
#include "/home/samir/Martini3/martini_v3.0.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_ions_v1.itp"
#include "/home/samir/Martini3/martini_v3.0.0_phospholipids_v1.itp"
#include "/home/samir/Martini3/martini_v3.0_sterols_v1.0.itp"
#include "/home/samir/Martini3/martini_v3.0.0_solvents_v1.itp"
EOF
cat topol.add sol.top > tmp
rm topol.add
mv tmp sol.top
sed -i -e 's/.*Insanely.*/Bicelle/g' sol.top
sed -i -e "s/.*Protein.*/POPC          $mol_POPC/g" sol.top
sed -i "/POPC/a\CHOL          $mol_CHOL" sol.top
sed -i "/POPC/a\POPS          $mol_POPS" sol.top
sed -i -e 's/NA+/NA /g' sol.top
sed -i -e 's/CL-/CL /g' sol.top

# change structure file
sed -i -e 's/NA+ /NA  /g' sol.gro
sed -i -e 's/ NA+/  NA/g' sol.gro
sed -i -e 's/CL- /CL  /g' sol.gro
sed -i -e 's/ CL-/  CL/g' sol.gro

# 2nd minimization
gmx grompp -f min.mdp -c sol.gro -p sol.top -o bc_min.tpr -quiet
gmx mdrun -deffnm bc_min -v -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet

# encurv
cat << EOF > plumed3.dat
# Define dummy atom at the center of curvature
piv: FIXEDATOM AT=35,35,12

# Atoms within 8nm of the dummy atom in the membrane
membr: GROUP NDX_FILE=bc_index.ndx NDX_GROUP=patch

# center of mass for patch
c1: COM ATOMS=membr

# Distance between centers of mass and virtual atom
d1: ZDISTANCES ATOMS1=piv,c1 LOWEST

# Define upper wall
uwal1: UPPER_WALLS ARG=d1.lowest AT=8 KAPPA=150

# Define lower wall
lwal1: UPPER_WALLS ARG=d1.lowest AT=6 KAPPA=150
EOF

# 2nd equilibration
cat << EOF > index.input
2 | 3 | 4
name 7 LIP
5 | 6
name 8 SOL
q
EOF
gmx make_ndx -f bc_min.gro -quiet -o bc_index.ndx < index.input
rm index.input
cat bc_index.ndx > tmp
mv tmp bc_index.ndx
rm select.ndx
gmx grompp -f bc_eq.mdp -c bc_min.gro -r bc_min.gro -p sol.top -o bc_eq.tpr -n bc_index.ndx -quiet
gmx mdrun -deffnm bc_eq -v -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# simulation
echo 0 | gmx trjconv -f bc_eq.xtc -s bc_eq.tpr -fr frame.ndx -o frame.gro -quiet
rm bc_index.ndx
cat << EOF > index.input
2 | 3 | 4
name 7 LIP
5 | 6
name 8 SOL
q
EOF
gmx make_ndx -f bc_min.gro -quiet -o bc_index.ndx < index.input
rm index.input
gmx select -f frame.gro -selrpos res_com -s bc_eq.tpr -quiet -on select.ndx -select 'same residue as (resname POPC or resname POPS or resname CHOL) and (within 10 of [33.48, 33.48, 10] or within 10 of [33.48, 33.48, 7.5] or within 10 of [33.48, 33.48, 12.5])'
sed -i -e 's/.*same_residue_as.*/[ patch ]/g' select.ndx
cat bc_index.ndx select.ndx > bc_index2.ndx
rm \select*
gmx grompp -f md.mdp -p sol.top -c frame.gro -n bc_index2.ndx -quiet -o md.tpr
gmx mdrun -v -deffnm md -plumed plumed3.dat -cpi md.cpt -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

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
fatslim apl -c eq.gro -t md.xtc -n index2.ndx --hg-group headgroups --export-apl-raw md_apl.csv
python3 show_apl_map.py
