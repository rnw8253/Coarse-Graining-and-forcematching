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
	global bond_fc,bond_equil_values,angle_fc,angle_equil_values,dihedral_fc,dihedral_period,dihedral_phase,nbonh,nbona,ntheta,ntheth,nphia,nphih,bondsh,bondsa,anglesh,anglesa,dihedralsh,dihedralsa,n_atoms,n_types,atom_names,atom_type_index,nb_parm_index,lj_a_coeff,lj_b_coeff,atom_types
	
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
			#dihedralsh = np.zeros((nphih,5),dtype=int)
			#dihedralsa = np.zeros((nphia,5),dtype=int)
			dihedralsh = []
			dihedralsa = []
			dihedralsh_temp = np.zeros((nphih,5),dtype=int)
			dihedralsa_temp = np.zeros((nphia,5),dtype=int)
			atom_names = []
			atom_types = []
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
				dihedralsh_temp[j][0] = dihedralsh_linear[j*5]
				dihedralsh_temp[j][1] = dihedralsh_linear[j*5+1]
				dihedralsh_temp[j][2] = dihedralsh_linear[j*5+2]
				dihedralsh_temp[j][3] = dihedralsh_linear[j*5+3]
				dihedralsh_temp[j][4] = dihedralsh_linear[j*5+4]
		
			for j in range(len(dihedralsh_temp)):
				if j == 0:
					dihedralsh.append([dihedralsh_temp[0][0],dihedralsh_temp[0][1],dihedralsh_temp[0][2],dihedralsh_temp[0][3],dihedralsh_temp[0][4]])
				else:
					if np.abs(dihedralsh_temp[j][0]) != np.abs(dihedralsh[-1][0]) or np.abs(dihedralsh_temp[j][1]) != np.abs(dihedralsh[-1][1]) or np.abs(dihedralsh_temp[j][2]) != np.abs(dihedralsh[-1][2]) or np.abs(dihedralsh_temp[j][3]) != np.abs(dihedralsh[-1][3]):
						dihedralsh.append([dihedralsh_temp[j][0],dihedralsh_temp[j][1],dihedralsh_temp[j][2],dihedralsh_temp[j][3],dihedralsh_temp[j][4]])
						
		if lines[i][0:32] == '%FLAG DIHEDRALS_WITHOUT_HYDROGEN':
			for j in range(n_dihedralsa_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					dihedralsa_linear[j*10+k] = int(temp[k])			
			for j in range(nphia):
				dihedralsa_temp[j][0] = dihedralsa_linear[j*5]
				dihedralsa_temp[j][1] = dihedralsa_linear[j*5+1]
				dihedralsa_temp[j][2] = dihedralsa_linear[j*5+2]
				dihedralsa_temp[j][3] = dihedralsa_linear[j*5+3]
				dihedralsa_temp[j][4] = dihedralsa_linear[j*5+4]
			
			for j in range(len(dihedralsa_temp)):
				if j == 0:
					dihedralsa.append([dihedralsa_temp[0][0],dihedralsa_temp[0][1],dihedralsa_temp[0][2],dihedralsa_temp[0][3],dihedralsa_temp[0][4]])
					
				else:
					if np.abs(dihedralsa_temp[j][0]) != np.abs(dihedralsa[-1][0]) or np.abs(dihedralsa_temp[j][1]) != np.abs(dihedralsa[-1][1]) or np.abs(dihedralsa_temp[j][2]) != np.abs(dihedralsa[-1][2]) or np.abs(dihedralsa_temp[j][3]) != np.abs(dihedralsa[-1][3]):
						dihedralsa.append([dihedralsa_temp[j][0],dihedralsa_temp[j][1],dihedralsa_temp[j][2],dihedralsa_temp[j][3],dihedralsa_temp[j][4]])
						
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

		if lines[i][0:21] == '%FLAG AMBER_ATOM_TYPE':
			for j in range(n_name_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					atom_types.append(temp[k])

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
out_file = sys.argv[2]

u = mda.Universe(top_file)
#u2 = mda.Universe(psf_file)
ParsePrmtopBonded(top_file)

################################
##### Prmtop to psf top #####
################################
#types_index_name = [[1,'GN'],[2,'GHN'],[3,'GCA'],[4,'GHC'],[5,'GCT'],[6,'OXT'],[7,'PHC4'],[8,''],[9,''],[10,''],[11,''],[12,'']]

out = open(out_file,'w')

# Write out the header
out.write("PSF\n\n")
out.write("       2 !NTITLE\n")
out.write(" REMARKS original generated structure x-plor psf file\n")
out.write(" REMARKS segment A { first NONE; last NONE; auto angles dihedrals }\n")

# Write out the atoms
out.write("\n   %s !NATOM\n" %(n_atoms))


for i in range(n_atoms):
	out.write("%8s A    %4s %4s %4s %4s %14.6f%14.6f%8s\n" %(i+1,u.atoms[i].resid,u.atoms[i].resname,u.atoms[i].name,atom_types[i],0,u.atoms[i].mass,0)) 

out.close

# Write out the bonds
n_bonds = nbona + nbonh

out.write("\n %7s !NBOND: bonds\n" %(n_bonds))
count = 1
for i in range(len(bondsh)):
	out.write(" %7s %7s" %(bondsh[i][0]/3 + 1,bondsh[i][1]/3 + 1))
	if count == 4:
		out.write("\n")
		count = 1
	else:
		count += 1

for i in range(len(bondsa)):
	out.write(" %7s %7s" %(bondsa[i][0]/3 + 1,bondsa[i][1]/3 + 1))
	if count == 4:
		out.write("\n")
		count = 1
	else:
		count += 1

# Write out the angles
n_angles = ntheta + ntheth

out.write("\n %7s !NTHETA: angles\n" %(n_angles))
count = 1
for i in range(len(anglesh)):
	out.write(" %7s %7s %7s" %(anglesh[i][0]/3 + 1,anglesh[i][1]/3 + 1,anglesh[i][2]/3 + 1))
	if count == 3:
		out.write("\n")
		count = 1
	else:
		count += 1

for i in range(len(anglesa)):
	out.write(" %7s %7s %7s" %(anglesa[i][0]/3 + 1,anglesa[i][1]/3 + 1,anglesa[i][2]/3 + 1))
	if count == 3:
		out.write("\n")
		count = 1
	else:
		count += 1

# Write out the dihedrals
n_dihedrals = len(dihedralsh) + len(dihedralsa)

out.write("\n %7s !NPHI: dihedrals\n" %(n_dihedrals))
count = 1
extra_count = 0
for i in range(len(dihedralsh)):
	out.write(" %7s %7s %7s %7s" %(np.abs(dihedralsh[i][0])/3 + 1,np.abs(dihedralsh[i][1])/3 + 1,np.abs(dihedralsh[i][2])/3 + 1,np.abs(dihedralsh[i][3])/3 + 1))
	if count == 2:
		out.write("\n")
		count = 1
	else:
		count += 1

for i in range(len(dihedralsa)):
	out.write(" %7s %7s %7s %7s" %(np.abs(dihedralsa[i][0])/3 + 1,np.abs(dihedralsa[i][1])/3 + 1,np.abs(dihedralsa[i][2])/3 + 1,np.abs(dihedralsa[i][3])/3 + 1))
	if count == 2:
		out.write("\n")
		count = 1
	else:
		count += 1

# Write out the imporopers
n_impropers = 0
n_donors = 0
n_acceptors = 0
n_nnbs = 0
n_ngrps1 = 1
n_ngrps2 = 0

out.write("\n %7s !NIMPHI: impropers\n" %(n_impropers))

out.write("\n\n %7s !NDON: donors\n" %(n_donors))

out.write("\n\n %7s !NACC: acceptors\n" %(n_acceptors))

out.write("\n\n %7s !NNB\n\n" %(n_nnbs))
for i in range(135):
	out.write("       0       0       0       0       0       0       0       0\n")

out.write("\n\n %7s %7s !NGRP:\n" %(n_ngrps1,n_ngrps2))
out.write("       0       0       0")

out.write("\n")
out.write("\n")

################################
