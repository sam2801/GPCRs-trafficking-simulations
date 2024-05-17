A repository for the bachelor thesis "An in silico study of how Helix 8 in G protein-coupled receptors affect trafficking through endocytosis".
Author: Samir Mustafa
17/05/2024

The Martini 3 folder includes the itp files that were used. 
The other folders contain a script running their particular simulation and their mdp files as well as any Python script that was used with it:
Bumpy_sim is the simulation that uses BUMPy to curve the bilayer.
EnCurv_2D is the simulation that uses EnCurv to induce a 2D curvature.
EnCurv_sim is the simulation that uses EnCurv to induce a 3D curvature.
bicelle_patch_center is the simulation that cuts off the edges of the bilayer and uses a Plumed distance wall on the center of the bilayer patch to prevent it from curving too much.
bicelle_patch_edges is the simulation that cuts off the edges of the bilayer and uses a Plumed distance wall on the centers of mass for areas around the edge to keep them and place and restrict the curvature of the central area.
insane_hole_gromacs is the simulation where an artificial pore in the center of the bilayer was made and maintained to alleviate the difference in lipid density between leaflets.
pepcurv is the simulation that embeds peptides into the membrane so they can induce curvature.
