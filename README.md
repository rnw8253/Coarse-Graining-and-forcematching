# Coarse-Graining-and-forcematching

makecg_30.tcl, make_cg_30_pdb.tcl and makecg_mono_pdb.tcl can all be used to cg the atoms in vmd from all-tom simulations and pdbs. These should be used if you are not using an all-atom cg where you are cg only the solvent. res_sel.py can be used to help make the selections from which you can copy and paste into the tcl scripts

using the output pdbs and trajectories, use top_from_pdb.py to create topology file and eventually a psf file

if you are just cg out the solvent then just use the prmtop_tp_psf.py script


Once you have a PSF file use imin_force_calc.sh and LJ2_q1.00.frc.nb.in to calculate the forces for an amber trajectory and out put an amber trajectory with forces
The settings in LJ2_q1.00.frc.nb.in can be changed depending on what forces you want. The current setting is for Lennard-jones forces only so that we can test the forcematching scripts and to make sure the bonds are defined correctly in the psf and also that we have enough sampling.

Using the force trajectory and psf file we can run the forcematch.py scripts to calculate nonbonded forces between atom-types

plot_lj.py and plot_lj_forces.py can be used to compare the calculed LJ potentials vs the ones defined in the prmtop 


