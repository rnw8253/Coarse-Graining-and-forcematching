import sys
import os
import numpy as np
import math
from MDAnalysis import Universe
from ForcePy import *
import pickle
#np.set_printoptions(threshold=sys.maxsize)

	
##############################
#######  Main Script  ########
##############################

# read in command line argument                                                                                                                                                                                                                                                                                          
top_file = sys.argv[1]
traj_file = sys.argv[2]

print "Topology file:", top_file
print "Trajectory file:", traj_file

# Declare universe
fine_uni = Universe(top_file, traj_file)
fine_uni.trajectory.periodic = True



#coarse_uni =CGUniverse(fine_uni,selections=['resname GC1 and name N','resname GC1 and name H','resname GC1 and name CA','resname GC1 and name HA2','resname GC1 and name HA3','resname GC1 and name C','resname GC1 and name O','resname GC1 and name OXT','resname G1 and name N','resname G1 and name H','resname G1 and name CA','resname G1 and name HA2','resname G1 and name HA3','resname G1 and name C','resname G1 and name O','resname PDI and name NA','resname PDI and name C8','resname PDI and name O1','resname PDI and name C9','resname PDI and name O2','resname PDI and name C10','resname PDI and name C11','resname PDI and name H16','resname PDI and name C12','resname PDI and name C13','resname PDI and name C14','resname PDI and name C15','resname PDI and name C16','resname PDI and name H17','resname PDI and name C17','resname PDI and name H18','resname PDI and name C18','resname PDI and name C19','resname PDI and name C20','resname PDI and name C21','resname PDI and name H19','resname PDI and name H20','resname PDI and name C22','resname PDI and name C23','resname PDI and name C24','resname PDI and name H21','resname PDI and name C25','resname PDI and name H22','resname PDI and name C26','resname PDI and name C27','resname PDI and name O3','resname PDI and name C28','resname PDI and name C29','resname PDI and name C30','resname PDI and name O4','resname PDI and name NB','resname PDI and name C31','resname PDI and name H23','resname PDI and name CAT','resname PDI and name OAT','resname PDI and name CA1','resname PDI and name HA1','resname PDI and name HA2','resname PDI and name CB1','resname PDI and name HB1','resname PDI and name HB2','resname PDI and name CBT','resname PDI and name OBT','resname G2 and name N','resname G2 and name H','resname G2 and name CA','resname G2 and name HA2','resname G2 and name HA3','resname G2 and name C','resname G2 and name O','resname GC2 and name N','resname GC2 and name H','resname GC2 and name CA','resname GC2 and name HA2','resname GC2 and name HA3','resname GC2 and name C','resname GC2 and name O','resname GC2 and name OXT'],names=['N','H','CA','HA2','HA3','CT','O','OXT','N','H','CA','HA2','HA3','C','O','NA','C8','O1','C9','O2','C10','C11','H16','C12','C13','C14','C15','C16','H17','C17','H18','C18','C19','C20','C21','H19','H20','C22','C23','C24','H21','C25','H22','C26','C27','O3','C28','C29','C30','O4','NB','C31','H23','CAT','OAT','CA1','HA1','HA2','CB1','HB1','HB2','CBT','OBT','N','H','CA','HA2','HA3','C','O','N','H','CA','HA2','HA3','CT','O','OXT'] ,types=['GN','GHN','GCA','GHC','GHC','GCT','OXT','OXT','GN','GHN','GCA','GHC','GHC','GC','GO','PN','PCO','PO','PCO','PO','PC1','PC4','PHC4','PC3','PC2','PC2','PC1','PC5','PHC5','PC4','PHC4','PC3','PC3','PC4','PC5','PHC5','PHC4','PC2','PC3','PC4','PHC4','PC5','PHC5','PC1','PCO','PO','PC2','PC1','PCO','PO','PN','PC5','PHC5','GC','GO','GCA','GHC','GHC','GCA','GHC','GHC','GC','GO','GN','GHN','GCA','GHC','GHC','GC','GO','GN','GHN','GCA','GHC','GHC','GCT','OXT','OXT'], collapse_hydrogens=False)
#
##coarse_uni = coarse_uni.cache()
#
#nres_mol = 9 # number of residues in a molecule
#
#add_residue_bonds(coarse_uni, 'name OXT', 'name CT')
#add_residue_bonds(coarse_uni, 'name O', 'name CT')
#add_residue_bonds(coarse_uni, 'name CT', 'name CA')
#add_residue_bonds(coarse_uni, 'name CA', 'name N')
#add_residue_bonds(coarse_uni, 'name CA', 'name HA2')
#add_residue_bonds(coarse_uni, 'name CA', 'name HA3')
#add_residue_bonds(coarse_uni, 'name N', 'name H')
#add_residue_bonds(coarse_uni, 'name C', 'name O')
#add_residue_bonds(coarse_uni, 'name C', 'name CA')
#
#add_residue_bonds(coarse_uni, 'name CAT', 'name OAT')
#add_residue_bonds(coarse_uni, 'name CAT', 'name CA1')
#add_residue_bonds(coarse_uni, 'name CA1', 'name HA1')
#add_residue_bonds(coarse_uni, 'name CA1', 'name HA2')
#add_residue_bonds(coarse_uni, 'name CA1', 'name NA')
#add_residue_bonds(coarse_uni, 'name NA', 'name C9')
#add_residue_bonds(coarse_uni, 'name NA', 'name C8')
#add_residue_bonds(coarse_uni, 'name C9', 'name O2')
#add_residue_bonds(coarse_uni, 'name C9', 'name C10')
#add_residue_bonds(coarse_uni, 'name C8', 'name O1')
#add_residue_bonds(coarse_uni, 'name C8', 'name C15')
#add_residue_bonds(coarse_uni, 'name C10', 'name C31')
#add_residue_bonds(coarse_uni, 'name C10', 'name C14')
#add_residue_bonds(coarse_uni, 'name C11', 'name H16')
#add_residue_bonds(coarse_uni, 'name C11', 'name C12')
#add_residue_bonds(coarse_uni, 'name C31', 'name C11')
#add_residue_bonds(coarse_uni, 'name C31', 'name H23')
#add_residue_bonds(coarse_uni, 'name C14', 'name C13')
#add_residue_bonds(coarse_uni, 'name C13', 'name C12')
#add_residue_bonds(coarse_uni, 'name C13', 'name C18')
#add_residue_bonds(coarse_uni, 'name C15', 'name C14')
#add_residue_bonds(coarse_uni, 'name C15', 'name C16')
#add_residue_bonds(coarse_uni, 'name C16', 'name C17')
#add_residue_bonds(coarse_uni, 'name C16', 'name H17')
#add_residue_bonds(coarse_uni, 'name C17', 'name C18')
#add_residue_bonds(coarse_uni, 'name C17', 'name H18')
#add_residue_bonds(coarse_uni, 'name C12', 'name C23')
#add_residue_bonds(coarse_uni, 'name C23', 'name C24')
#add_residue_bonds(coarse_uni, 'name C23', 'name C22')
#add_residue_bonds(coarse_uni, 'name C18', 'name C19')
#add_residue_bonds(coarse_uni, 'name C19', 'name C22')
#add_residue_bonds(coarse_uni, 'name C19', 'name C20')
#add_residue_bonds(coarse_uni, 'name C24', 'name C25')
#add_residue_bonds(coarse_uni, 'name C24', 'name H21')
#add_residue_bonds(coarse_uni, 'name C25', 'name C26')
#add_residue_bonds(coarse_uni, 'name C25', 'name H22')
#add_residue_bonds(coarse_uni, 'name C26', 'name C27')
#add_residue_bonds(coarse_uni, 'name C22', 'name C28')
#add_residue_bonds(coarse_uni, 'name C28', 'name C26')
#add_residue_bonds(coarse_uni, 'name C28', 'name C29')
#add_residue_bonds(coarse_uni, 'name C20', 'name H20')
#add_residue_bonds(coarse_uni, 'name C20', 'name C21')
#add_residue_bonds(coarse_uni, 'name C21', 'name C29')
#add_residue_bonds(coarse_uni, 'name C21', 'name H19')
#add_residue_bonds(coarse_uni, 'name C29', 'name C30')
#add_residue_bonds(coarse_uni, 'name C27', 'name O3')
#add_residue_bonds(coarse_uni, 'name C27', 'name NB')
#add_residue_bonds(coarse_uni, 'name C30', 'name O4')
#add_residue_bonds(coarse_uni, 'name C30', 'name NB')
#add_residue_bonds(coarse_uni, 'name NB', 'name CB1')
#add_residue_bonds(coarse_uni, 'name CB1', 'name HB1')
#add_residue_bonds(coarse_uni, 'name CB1', 'name HB2')
#add_residue_bonds(coarse_uni, 'name CB1', 'name CBT')
#add_residue_bonds(coarse_uni, 'name CBT', 'name OBT')
#
#
#add_sequential_bonds(coarse_uni,nres_mol) # function that adds a bond between every residue in a molecule 

write_structure(fine_uni, 'cg_1.pdb', bonds='all')

fm = ForceMatch(fine_uni) 


ff = FileForce() #This just says the forces are found in the universe passed to the ForceMatch object
fm.add_ref_force(ff)  #Set this force as the reference force to be matched  

pair_mesh = Mesh.UniformMesh(0,12,0.05) #This is the mesh on which the force-field will be built. It is in Angstroms
pairwise_force = SpectralForce(Pairwise, pair_mesh, Basis.UnitStep)  

fm.add_and_type_pair(pairwise_force) #Copy this force type and clone it for each pair-interaction type


fm.force_match()
#fm.force_match_mpi()

fm.write_lammps_scripts()





