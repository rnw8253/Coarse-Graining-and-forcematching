#USAGE :  python python_mdanalysis_skeleton.py [config file name]

#DEPENDCIES : numpy, MDAnalysis, math

#CONFIG FILE FORMAT:
#   TopFile = [topology file name (prmtop file)]
#   TrajFile = [trajectory file name (mdcrd file)]
#   OutFile = [output data file name]

import sys
import os
import numpy as np
import math
import MDAnalysis

# read the configuration file and populate the global variables
def ParsePdbFile(cgpdb_file):
	global atom_name, res_name, chain_name, res_num, atom_coord, n_atoms, x_atom_name
	f = open(cgpdb_file)
	atom_name = []
	res_name = []
	chain_name = []
	res_num = []
	atom_coord = []
	n_atoms = 0
	for line in f:
		line_list=line.split()
		if line_list[0] == 'ATOM':
			atom_name.append(line_list[2])
			res_name.append(line_list[3])
			chain_name.append(line_list[4])
			res_num.append(line_list[5])
			atom_coord.append([])
			atom_coord[n_atoms].append(line_list[6])
			atom_coord[n_atoms].append(line_list[7])
			atom_coord[n_atoms].append(line_list[8])
			n_atoms += 1
			
		# first remove comments
		#if '#' in line:
		#	line, comment = line.split('#',1)
		#if 'ATOM' in line:
			#atom_name.append(line[13:16])
			#res_name.append(line[17:20])
			#chain_name.append(line[21:22])
			#res_num.append(line[23:26])
			#atom_coord.append([])
			#atom_coord[n_atoms].append(line[30:38])
			#atom_coord[n_atoms].append(line[38:46])
			#atom_coord[n_atoms].append(line[46:54])
			#n_atoms += 1
	f.close()
	atom_coord = np.asmatrix(atom_coord,dtype=float)
	x_atom_name = []
	for i in range(len(atom_name)):
		x_atom_name.append('+%s' %(atom_name[i]))
# create psf topology file
def CreateTopFile(top_file):
	global atom_name, res_name, chain_name, res_num, atom_coord, n_atoms
	f = open(top_file, 'w')

	atom_type_index_count = 1
	for i in range(n_CG_particles):
		if atom_type_index[i] == atom_type_index_count:
			mass = 0
			for j in range(len(CG_par_sel[i].atoms)):
				if CG_par_sel[i].atoms[j].mass == 22.98977:
					mass += 14.007
				else:
					mass += CG_par_sel[i].atoms[j].mass
			f.write("MASS     %2s  %4s %7.3f\n" %(atom_type_index_count,atom_type[i],mass)) 
			atom_type_index_count += 1
			
	f.write("\n\nAUTO ANGLES DIHE\n\n")
	res_type_count = 0
	resnum_count = 0
	count = 0
	count1 = 0
	for j in range(n_CG_particles):
#		print j
		if int(res_num[j]) != resnum_count:
			if res_type_index[j] - res_type_count == 1:
				if j != 0:
					print count1, j-1
					for i in range(count1,j):
						if bond_index[i][0] != 0:
							if i == count1:
								f.write("BOND %4s %4s" %(atom_name[i], atom_name[bond_index[i][0]]))
							else:
								f.write("\nBOND %4s %4s" %(atom_name[i], atom_name[bond_index[i][0]]))
							if bond_index[i][1] != 0:
								f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][1]]))
								if bond_index[i][2] != 0:
									f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][2]]))
						if bond_index[i][3] != 0:
							f.write("   %4s %4s" %(atom_name[i], x_atom_name[bond_index[i][3]]))
					f.flush()
					count1 = j
	
				f.write("\n\n")
				f.write("RESI %3s       0.000\n" %(res_name[j]))
				f.write("GROUP\n")
				f.write("ATOM %4s %4s    %6.3f\n" % (atom_name[j], atom_type[j], 0.0))
				res_type_count += 1
				resnum_count = count + 1
			count = int(res_num[j]) 
			
		else:
			if res_type_index[j] - res_type_count == 0:
					f.write("ATOM %4s %4s    %6.3f\n" % (atom_name[j], atom_type[j], 0.0))

		if j == 39:
			for i in range(count1,j):
				if bond_index[i][0] != 0:
					if i == count1:
						f.write("BOND %4s %4s" %(atom_name[i], atom_name[bond_index[i][0]]))
					else:
						f.write("\nBOND %4s %4s" %(atom_name[i], atom_name[bond_index[i][0]]))
					if bond_index[i][1] != 0:
						f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][1]]))
						if bond_index[i][2] != 0:
							f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][2]]))
				if bond_index[i][3] != 0:
					f.write("   %4s %4s" %(atom_name[i], x_atom_name[bond_index[i][3]]))
		
#	f.write("\n\n")

#	for i in range(n_CG_particles):
#		num_bonds = len(bond_index[i])
#		if bond_index[i][0] != 0:
#			f.write("BOND %4s %4s" %(atom_name[i], atom_name[bond_index[i][0]]))
#			if bond_index[i][1] != 0:
#				f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][1]]))
#				if bond_index[i][2] != 0:
#					f.write("   %4s %4s" %(atom_name[i], atom_name[bond_index[i][2]]))
#		if bond_index[i][3] != 0:
#			f.write("   %4s +%s\n" %(atom_name[i], atom_name[bond_index[i][3]]))
#		else:
#			f.write("\n")
#				
#		#	f.write("BOND %4s %4s   %4s %4s   %4s %4s\n" % (atom_name[bp], atom_name[bp+1], atom_name[bp], atom_name[n_atoms-bp-1], atom_name[n_atoms-bp-1], atom_name[n_atoms-bp-2]))
#		f.flush()
	
	#f.write("BOND %4s %4s\n\n" %(atom_name[n_bps-1], atom_name[n_bps]))
	f.close()

def CreatePsfFile(psf_file):
	global res_name, pdb_file, top_file
	
	pgn_file = "GGGGG.pgn"
	new_pdb_file = "GGGGG.cg.v02.pdb"

	f = open(pgn_file, 'w')
	f.write("topology %s\n" % (top_file))
	f.write("segment A {\n")
	f.write("   first NONE\n")
	f.write("   last NONE\n")
	f.write("   pdb %s\n" % (cgpdb_file))
	f.write("}\n")
	f.write("coordpdb %s   A\n" % (cgpdb_file))
	f.write("writepsf %s\n" % (psf_file))
	f.write("writepdb %s\n" % (new_pdb_file))

	f.close()

	command_string = "psfgen " + pgn_file
	os.system(command_string)

# read in command line argument
cgpdb_file = sys.argv[1]
#psf_file = sys.argv[2]
pdb_file = sys.argv[2]
u = MDAnalysis.Universe(pdb_file)
CG_particle_list = ["resid 1 and name OXT","resid 1 and name O","resid 1 and name C","resid 1 and name CA","resid 1 and name HA2","resid 1 and name HA3","resid 1 and name N","resid 1 and name H","resid 2 and name C","resid 2 and name O","resid 2 and name CA","resid 2 and name HA2","resid 2 and name HA3","resid 2 and name N","resid 2 and name H","resid 3 and name C","resid 3 and name O","resid 3 and name CA","resid 3 and name HA2","resid 3 and name HA3","resid 3 and name N","resid 3 and name H","resid 4 and name C","resid 4 and name O","resid 4 and name CA","resid 4 and name HA2","resid 4 and name HA3","resid 4 and name N","resid 4 and name H","resid 5 and name CA","resid 5 and name OA","resid 5 and name CA1","resid 5 and name HA1","resid 5 and name HA2","resid 5 and name NA","resid 5 and name C9","resid 5 and name O2","resid 5 and name C8","resid 5 and name O1","resid 5 and name C10","resid 5 and name C11","resid 5 and name C31","resid 5 and name H23","resid 5 and name H16","resid 5 and name C14","resid 5 and name C13","resid 5 and name C15","resid 5 and name C16","resid 5 and name C17","resid 5 and name H17","resid 5 and name H18","resid 5 and name C12","resid 5 and name C23","resid 5 and name C18","resid 5 and name C19","resid 5 and name C24","resid 5 and name C25","resid 5 and name C26","resid 5 and name H21","resid 5 and name H22","resid 5 and name C22","resid 5 and name C28","resid 5 and name C20","resid 5 and name C21","resid 5 and name C29","resid 5 and name H19","resid 5 and name H20","resid 5 and name C27","resid 5 and name O3","resid 5 and name C30","resid 5 and name O4","resid 5 and name NB","resid 5 and name CB1","resid 5 and name HB1","resid 5 and name HB2","resid 5 and name CB","resid 5 and name OB","resid 6 and name N","resid 6 and name H","resid 6 and name CA","resid 6 and name HA2","resid 6 and name HA3","resid 6 and name C","resid 6 and name O","resid 7 and name N","resid 7 and name H","resid 7 and name CA","resid 7 and name HA2","resid 7 and name HA3","resid 7 and name C","resid 7 and name O","resid 8 and name N","resid 8 and name H","resid 8 and name CA","resid 8 and name HA2","resid 8 and name HA3","resid 8 and name C","resid 8 and name O","resid 9 and name N","resid 9 and name H","resid 9 and name CA","resid 9 and name HA2","resid 9 and name HA3","resid 9 and name OXT","resid 9 and name O","resid 9 and name C"]

n_CG_particles = len(CG_particle_list)
CG_par_sel = []
for i in range(n_CG_particles):
	CG_par_sel.append(u.select_atoms(CG_particle_list[i]))

#              0     1      2     3     4      5     6      7     8     9     10    11    12   13    14    15    16     17    18    19    20    21    22   23     24    25    26    27    28    29     30    31    32     33    34    35     36    37    38    39
atom_type = ["OXT","O","C","C","H","H","N","H"]
atom_type_index = [1,2,3,4,5,6,4,5,6,4,5,6,7,8,9,9,10,11,10,12,12,10,11,10,9,9,8,7,6,5,4,6,5,4,6,5,4,3,2,1] # 1 indexed
res_type_index = [1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5] # 1 indexed
bond_index=[[1,0,0,0],[2,0,0,0],[0,0,0,3],[4,0,0,0],[5,0,0,0],[0,0,0,6],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,12],[13,0,0,0],[14,15,0,0],[16,0,0,0],[18,0,0,0],[17,19,0,0],[18,19,20,0],[20,0,0,0],[21,22,0,0],[22,23,0,0],[22,24,0,0],[23,0,0,0],[25,0,0,0],[26,0,0,0],[26,0,0,0],[27,0,0,0],[0,0,0,28],[29,0,0,0],[30,0,0,0],[0,0,0,31],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[38,0,0,0],[39,0,0,0],[0,0,0,0]]
# read pdb file
ParsePdbFile(cgpdb_file)
top_file = "GGGGG.cg.top"
# create topology file
CreateTopFile(top_file)

# create pgn file
CreatePsfFile('GGGGG.cg.psf')

# 
