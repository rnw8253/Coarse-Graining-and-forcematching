import sys
import os
import numpy as np
import math
import MDAnalysis as mda
from MDAnalysis.analysis.align import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.ticker import NullFormatter

def ParsePrmtopBonded(top_file):
	global bond_fc,bond_equil_values,angle_fc,angle_equil_values,dihedral_fc,dihedral_period,dihedral_phase,nbonh,nbona,ntheta,ntheth,nphia,nphih,bondsh,bondsa,anglesh,anglesa,dihedralsh,dihedralsa,n_atoms,n_types,atom_names,atom_type_index,nb_parm_index,lj_a_coeff,lj_b_coeff
	
	param = open(top_file,'r')
	pointers = np.zeros(31,dtype=np.int)
	lines = param.readlines()
	for i in range(len(lines)):	
		if lines[i][0:14] == '%FLAG POINTERS':
			for j in range(4):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					pointers[j*10+k] = int(temp[k])
			n_atoms = pointers[0]
			n_types = pointers[1]
			nbonh = pointers[2]
			nbona = pointers[12]
			ntheth = pointers[4]
			ntheta = pointers[13]
			nphih = pointers[6]
			nphia = pointers[14]
			numbnd = pointers[15]
			numang = pointers[16]
			numtra = pointers[17]
			n_type_lines = int(math.ceil(n_atoms/10.))
			n_name_lines = int(math.ceil(n_atoms/20.))
			n_nb_parm_lines = int(math.ceil(n_types*n_types/10.))
			n_lj_param_lines = int(math.ceil((n_types*(n_types+1)/2)/5.))
			n_bond_lines = int(math.ceil(numbnd/5.))
			n_angle_lines = int(math.ceil(numang/5.))
			n_dihedral_lines = int(math.ceil(numtra/5.))
			n_bondsh_lines = int(math.ceil(nbonh*3/10.))
			n_bondsa_lines = int(math.ceil(nbona*3/10.))
			n_anglesh_lines = int(math.ceil(ntheth*4/10.))
			n_anglesa_lines = int(math.ceil(ntheta*4/10.))
			n_dihedralsh_lines = int(math.ceil(nphih*5/10.))
			n_dihedralsa_lines = int(math.ceil(nphia*5/10.))
			bond_fc = np.zeros(numbnd,dtype=np.float)
			bond_equil_values = np.zeros(numbnd,dtype=np.float)
			angle_fc = np.zeros(numang,dtype=np.float)
			angle_equil_values = np.zeros(numang,dtype=np.float)
			dihedral_fc = np.zeros(numtra,dtype=np.float)
			dihedral_period = np.zeros(numtra,dtype=np.float)
			dihedral_phase = np.zeros(numtra,dtype=np.float)
			SCEE_factor = np.zeros(numtra,dtype=np.float)
			SCNB_factor = np.zeros(numtra,dtype=np.float)
			bondsh_linear = np.zeros(3*nbonh,dtype=int)
			bondsa_linear = np.zeros(3*nbona,dtype=int)
			bondsh = np.zeros((nbonh,3),dtype=int)
			bondsa = np.zeros((nbona,3),dtype=int)
			anglesh_linear = np.zeros(4*ntheth,dtype=int)
			anglesa_linear = np.zeros(4*ntheta,dtype=int)
			anglesh = np.zeros((ntheth,4),dtype=int)
			anglesa = np.zeros((ntheta,4),dtype=int)
			dihedralsh_linear = np.zeros(5*nphih,dtype=int)
			dihedralsa_linear = np.zeros(5*nphia,dtype=int)
			dihedralsh = np.zeros((nphih,5),dtype=int)
			dihedralsa = np.zeros((nphia,5),dtype=int)
			atom_names = []
			atom_type_index = np.zeros((n_atoms),dtype=int)
			nb_parm_index = np.zeros(n_types*n_types,dtype=int)
			lj_a_coeff = np.zeros((n_types*(n_types+1))/2,dtype=float)
			lj_b_coeff = np.zeros((n_types*(n_types+1))/2,dtype=float)


		if lines[i][0:25] == '%FLAG BOND_FORCE_CONSTANT':
			for j in range(n_bond_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					bond_fc[j*5+k] = float(temp[k])
		if lines[i][0:22] == '%FLAG BOND_EQUIL_VALUE':
			for j in range(n_bond_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					bond_equil_values[j*5+k] = float(temp[k])
		if lines[i][0:26] == '%FLAG ANGLE_FORCE_CONSTANT':
			for j in range(n_angle_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					angle_fc[j*5+k] = float(temp[k])
		if lines[i][0:23] == '%FLAG ANGLE_EQUIL_VALUE':
			for j in range(n_angle_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					angle_equil_values[j*5+k] = float(temp[k])
		if lines[i][0:29] == '%FLAG DIHEDRAL_FORCE_CONSTANT':
			for j in range(n_dihedral_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedral_fc[j*5+k] = float(temp[k])
		if lines[i][0:26] == '%FLAG DIHEDRAL_PERIODICITY':
			for j in range(n_dihedral_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedral_period[j*5+k] = float(temp[k])
		if lines[i][0:20] == '%FLAG DIHEDRAL_PHASE':
			for j in range(n_dihedral_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedral_phase[j*5+k] = float(temp[k])
		if lines[i][0:23] == '%FLAG SCEE_SCALE_FACTOR':
			for j in range(n_dihedral_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					SCEE_factor[j*5+k] = float(temp[k])
		if lines[i][0:23] == '%FLAG SCNB_SCALE_FACTOR':
			for j in range(n_dihedral_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					SCNB_factor[j*5+k] = float(temp[k])
		if lines[i][0:24] == '%FLAG BONDS_INC_HYDROGEN':
			for j in range(n_bondsh_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					bondsh_linear[j*10+k] = int(temp[k])
			for j in range(nbonh):
				bondsh[j][0] = bondsh_linear[j*3]
				bondsh[j][1] = bondsh_linear[j*3+1]
				bondsh[j][2] = bondsh_linear[j*3+2]
		if lines[i][0:28] == '%FLAG BONDS_WITHOUT_HYDROGEN':
			for j in range(n_bondsa_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					bondsa_linear[j*10+k] = int(temp[k])			
			for j in range(nbona):
				bondsa[j][0] = bondsa_linear[j*3]
				bondsa[j][1] = bondsa_linear[j*3+1]
				bondsa[j][2] = bondsa_linear[j*3+2]
		if lines[i][0:25] == '%FLAG ANGLES_INC_HYDROGEN':
			for j in range(n_anglesh_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					anglesh_linear[j*10+k] = int(temp[k])
			for j in range(ntheth):
				anglesh[j][0] = anglesh_linear[j*4]
				anglesh[j][1] = anglesh_linear[j*4+1]
				anglesh[j][2] = anglesh_linear[j*4+2]
				anglesh[j][3] = anglesh_linear[j*4+3]
		if lines[i][0:29] == '%FLAG ANGLES_WITHOUT_HYDROGEN':
			for j in range(n_anglesa_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					anglesa_linear[j*10+k] = int(temp[k])			
			for j in range(ntheta):
				anglesa[j][0] = anglesa_linear[j*4]
				anglesa[j][1] = anglesa_linear[j*4+1]
				anglesa[j][2] = anglesa_linear[j*4+2]
				anglesa[j][3] = anglesa_linear[j*4+3]
		if lines[i][0:28] == '%FLAG DIHEDRALS_INC_HYDROGEN':
			for j in range(n_dihedralsh_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedralsh_linear[j*10+k] = int(temp[k])
			for j in range(nphih):
				dihedralsh[j][0] = dihedralsh_linear[j*5]
				dihedralsh[j][1] = dihedralsh_linear[j*5+1]
				dihedralsh[j][2] = dihedralsh_linear[j*5+2]
				dihedralsh[j][3] = dihedralsh_linear[j*5+3]
				dihedralsh[j][4] = dihedralsh_linear[j*5+4]
		if lines[i][0:32] == '%FLAG DIHEDRALS_WITHOUT_HYDROGEN':
			for j in range(n_dihedralsa_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedralsa_linear[j*10+k] = int(temp[k])			
			for j in range(nphia):
				dihedralsa[j][0] = dihedralsa_linear[j*5]
				dihedralsa[j][1] = dihedralsa_linear[j*5+1]
				dihedralsa[j][2] = dihedralsa_linear[j*5+2]
				dihedralsa[j][3] = dihedralsa_linear[j*5+3]
				dihedralsa[j][4] = dihedralsa_linear[j*5+4]
		
		if lines[i][0:15] == '%FLAG ATOM_NAME':
			for j in range(n_name_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					atom_names.append(temp[k])

		if lines[i][0:21] == '%FLAG ATOM_TYPE_INDEX':
			for j in range(n_type_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					atom_type_index[j*10+k] = float(temp[k])

		if lines[i][0:26] == '%FLAG NONBONDED_PARM_INDEX':
			for j in range(n_nb_parm_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					nb_parm_index[j*10+k] = float(temp[k])


		if lines[i][0:25] == '%FLAG LENNARD_JONES_ACOEF':
			for j in range(n_lj_param_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					lj_a_coeff[j*5+k] = float(temp[k])
		if lines[i][0:25] == '%FLAG LENNARD_JONES_BCOEF':
			for j in range(n_lj_param_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					lj_b_coeff[j*5+k] = float(temp[k])

top_file = sys.argv[1]
pdb_file = sys.argv[2]
psf_file = sys.argv[3]

u = mda.Universe(top_file,pdb_file)
u2 = mda.Universe(psf_file)
ParsePrmtopBonded(top_file)

################################
##### Prmtop to lammps top #####
################################

out = open("GGGGG.30.cg.lammpstop",'w')

# Write out the number of atoms, bonds, angles, dihedrals, impropers, atom types, bond types, angle types, dihedrals types, improper types and box size
atom_types = []
for i in range(len(u2.atoms)):
	if i == 0:
		atom_types.append((u2.atoms[i].type,u2.atoms[i].mass)) 
	else:
		break_flag = False
		for j in range(len(atom_types)):
			if u2.atoms[i].type == atom_types[j][0]:
				break_flag = True
				break
		if break_flag == False:	
			atom_types.append((u2.atoms[i].type,u2.atoms[i].mass))


n_atom_types = len(atom_types)
n_bonds = len(u2.bonds)
n_bond_types = len(u2.bonds.types())
n_angles = len(u2.angles)
n_angle_types = len(u2.angles.types())
n_dihs = len(u2.dihedrals)
n_dih_types = len(u2.dihedrals.types())
#n_impropers = 0
n_impropers = len(u2.impropers)
if n_impropers == 0:
	n_improper_types = 0
else:
	n_improper_types = len(u2.impropers.types())

box = [120,120,120]
x_half = box[0]/2
y_half = box[1]/2
z_half = box[2]/2

out.write('    %4s atoms\n' %(n_atoms))
out.write('    %4s bonds\n' %(n_bonds))
out.write('    %4s angles\n' %(n_angles))
out.write('    %4s dihedrals\n' %(n_dihs))
out.write('    %4s impropers\n' %(n_impropers))
out.write('    %4s atom types\n' %(n_atom_types))
out.write('    %4s bond types\n' %(n_bond_types))
out.write('    %4s angle types\n' %(n_angle_types))
out.write('    %4s dihedral types\n' %(n_dih_types))
out.write('    %4s improper types\n' %(n_improper_types))
out.write('-%8.6f %8.6f xlo xhi\n' %(x_half,x_half))
out.write('-%8.6f %8.6f ylo yhi\n' %(y_half,y_half))
out.write('-%8.6f %8.6f zlo zhi\n' %(z_half,z_half))


# Write out the masses of each atom type
out.write('\n\n Masses\n\n')

for i in range(n_atom_types):
	out.write('       %2s  %6.3f # %4s\n' %(i+1,atom_types[i][1],atom_types[i][0]))

# Write out the atoms [count, res_num, atom_type_index, charge?, x-pos, y-pos, z-pos, #, name,  res_name]
out.write('\n Atoms \n\n')
for i in range(len(u2.atoms)):
	for j in range(len(atom_types)):
		if atom_types[j][0] == u2.atoms[i].type:
			out.write('%4s %3s %2s 0.000000 %9.6f %9.6f %9.6f # %s\n' %(i+1, u2.atoms[i].resid, j+1, u.atoms[i].position[0], u.atoms[i].position[1], u.atoms[i].position[2], u2.atoms[i].name))    
			break

# Write bonds [ count, type, atom1, atom2]
out.write('\n Bonds\n\n')

for i in range(len(u2.bonds)):
	for j in range(len(u2.bonds.types())):
		if u2.bonds[i].type == u2.bonds.types()[j]:
			out.write('%s %s %s %s\n' %(i+1, j+1, u2.bonds[i][0].index+1, u2.bonds[i][1].index+1))
			break

# Write angles [ count, type, atom1, atom2, atom3] 
out.write('\n Angles\n\n')

for i in range(len(u2.angles)):
	for j in range(len(u2.angles.types())):
		if u2.angles[i].type == u2.angles.types()[j]:
			out.write('%s %s %s %s %s\n' %(i+1, j+1, u2.angles[i][0].index+1, u2.angles[i][1].index+1, u2.angles[i][2].index+1))
			break
# Write dihedral [ count, type, atom1, atom2, atom3, atom4]
out.write('\n Dihedrals\n\n')

for i in range(len(u2.dihedrals)):
	for j in range(len(u2.dihedrals.types())):
		if u2.dihedrals[i].type == u2.dihedrals.types()[j]:
			out.write('%s %s %s %s %s %s\n' %(i+1, j+1, u2.dihedrals[i][0].index+1, u2.dihedrals[i][1].index+1, u2.dihedrals[i][2].index+1, u2.dihedrals[i][3].index+1))
			break


################################

################################
##### Create lammps par ########
################################

out2 = open("GGGGG.30.cg.lammpspar",'w')
# Write bond_coeff [ bond_coeff, bond_type, K (including 1/2 factor), r0 ]

for i in range(len(u2.bonds.types())):
	for j in range(len(u2.bonds)):
		if u2.bonds[j].type == u2.bonds.types()[i]:
			# find the atom indexs
			atoms1_index = u2.bonds[j][0].index + 1
			atoms2_index = u2.bonds[j][1].index + 1
			if u2.bonds[j][0].name[0] == 'H' or u2.bonds[j][1].name[0] == 'H':
				for k in range(len(bondsh)):
					if (bondsh[k][0] == (atoms1_index-1)*3 and bondsh[k][1] == (atoms2_index-1)*3) or (bondsh[k][0] == (atoms2_index-1)*3 and bondsh[k][1] == (atoms1_index-1)*3):
						constant_index = bondsh[k][2] - 1 
						out2.write('bond_coeff %2s %6.4f %6.4f                   # [%s] -- [%s]\n' %(i+1, float(bond_fc[constant_index]),  float(bond_equil_values[constant_index]), u2.bonds.types()[i][0], u2.bonds.types()[i][1]))
						break
			else:
				for k in range(len(bondsa)):
					if (bondsa[k][0] == (atoms1_index-1)*3 and bondsa[k][1] == (atoms2_index-1)*3) or (bondsa[k][0] == (atoms2_index-1)*3 and bondsa[k][1] == (atoms1_index-1)*3):
						constant_index = bondsa[k][2] - 1
						out2.write('bond_coeff %2s %6.4f %6.4f                   # [%s] -- [%s]\n' %(i+1, float(bond_fc[constant_index]),  float(bond_equil_values[constant_index]), u2.bonds.types()[i][0], u2.bonds.types()[i][1]))
						break
			break
		
# Write angle_coeff [ angle_coeff, angle_type, K (including 1/2 factor), theta0 ]

for i in range(len(u2.angles.types())):
	for j in range(len(u2.angles)):
		if u2.angles[j].type == u2.angles.types()[i]:
			# find the atom indexs
			atoms1_index = u2.angles[j][0].index + 1
			atoms2_index = u2.angles[j][1].index + 1
			atoms3_index = u2.angles[j][2].index + 1

			if u2.angles[j][0].name[0] == 'H' or u2.angles[j][1].name[0] == 'H' or u2.angles[j][2].name[0] == 'H':
				for k in range(len(anglesh)):
					if (anglesh[k][0] == (atoms1_index-1)*3 and anglesh[k][1] == (atoms2_index-1)*3 and anglesh[k][2] == (atoms3_index-1)*3) or (anglesh[k][0] == (atoms3_index-1)*3 and anglesh[k][1] == (atoms2_index-1)*3 and anglesh[k][2] == (atoms1_index-1)*3):
						constant_index = anglesh[k][3] - 1 
						out2.write('angle_coeff %2s %6.4f %6.4f                   # [%s] -- [%s] -- [%s]\n' %(i+1, float(angle_fc[constant_index]),  float(angle_equil_values[constant_index]), u2.angles.types()[i][0], u2.angles.types()[i][1], u2.angles.types()[i][2]))
						break
			else:
				for k in range(len(anglesa)):
					if (anglesa[k][0] == (atoms1_index-1)*3 and anglesa[k][1] == (atoms2_index-1)*3 and anglesa[k][2] == (atoms3_index-1)*3) or (anglesa[k][0] == (atoms3_index-1)*3 and anglesa[k][1] == (atoms2_index-1)*3 and anglesa[k][2] == (atoms1_index-1)*3):
						constant_index = anglesa[k][3] - 1
						out2.write('angle_coeff %2s %6.4f %6.4f                   # [%s] -- [%s] -- [%s]\n' %(i+1, float(angle_fc[constant_index]),  float(angle_equil_values[constant_index]), u2.angles.types()[i][0], u2.angles.types()[i][1], u2.angles.types()[i][2]))
						break
			break


# Write dihedral_coeff [ dihedral_coeff, dihedral_type, K, d (+1 or -1) , n (integer > 0)]

for i in range(len(u2.dihedrals.types())):
	for j in range(len(u2.dihedrals)):
		if u2.dihedrals[j].type == u2.dihedrals.types()[i]:
			# find the atom indexs
			atoms1_index = u2.dihedrals[j][0].index + 1
			atoms2_index = u2.dihedrals[j][1].index + 1
			atoms3_index = u2.dihedrals[j][2].index + 1
			atoms4_index = u2.dihedrals[j][3].index + 1

			if u2.dihedrals[j][0].name[0] == 'H' or u2.dihedrals[j][1].name[0] == 'H' or u2.dihedrals[j][2].name[0] == 'H' or u2.dihedrals[j][3].name[0] == 'H':
				for k in range(len(dihedralsh)):
					if (dihedralsh[k][0] == (atoms1_index-1)*3 and dihedralsh[k][1] == (atoms2_index-1)*3 and dihedralsh[k][2] == (atoms3_index-1)*3 and dihedralsh[k][3] == (atoms4_index-1)*3) or (dihedralsh[k][0] == (atoms4_index-1)*3 and dihedralsh[k][1] == (atoms3_index-1)*3 and dihedralsh[k][2] == (atoms2_index-1)*3 and dihedralsh[k][3] == (atoms1_index-1)*3):
						constant_index = dihedralsh[k][4] - 1 
						out2.write('dihedral_coeff %2s %6.4f %4s %6.4f            # [%s] -- [%s] -- [%s] -- [%s]\n' %(i+1, float(dihedral_fc[constant_index]), np.floor(np.cos(dihedral_phase[constant_index])),  float(dihedral_period[constant_index]), u2.dihedrals.types()[i][0], u2.dihedrals.types()[i][1], u2.dihedrals.types()[i][2], u2.dihedrals.types()[i][3]))
						break
			else:
				for k in range(len(dihedralsa)):
					if (dihedralsa[k][0] == (atoms1_index-1)*3 and dihedralsa[k][1] == (atoms2_index-1)*3 and dihedralsa[k][2] == (atoms3_index-1)*3 and dihedralsa[k][3] == (atoms4_index-1)*3) or (dihedralsa[k][0] == (atoms4_index-1)*3 and dihedralsa[k][1] == (atoms3_index-1)*3 and dihedralsa[k][2] == (atoms2_index-1)*3 and dihedralsa[k][3] == (atoms1_index-1)*3):
						constant_index = dihedralsa[k][4] - 1
						out2.write('dihedral_coeff %2s %6.4f %4s %6.4f            # [%s] -- [%s] -- [%s] -- [%s]\n' %(i+1, float(dihedral_fc[constant_index]), np.floor(np.cos(dihedral_phase[constant_index])), float(dihedral_period[constant_index]), u2.dihedrals.types()[i][0], u2.dihedrals.types()[i][1], u2.dihedrals.types()[i][2], u2.dihedrals.types()[i][3]))
						break
			break

################################
