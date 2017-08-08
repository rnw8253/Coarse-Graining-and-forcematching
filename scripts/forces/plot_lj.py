import sys
import os
import numpy as np
import math
import MDAnalysis
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
			print lj_a_coeff
		if lines[i][0:25] == '%FLAG LENNARD_JONES_BCOEF':
			for j in range(n_lj_param_lines):
				temp = lines[i+2+j].split()
				for k in range(len(temp)):
					lj_b_coeff[j*5+k] = float(temp[k])



def plot_1d(xdata1, ydata1, lj_a_coeff, lj_b_coeff, color1, color2, name1, name2):

        plt.plot(xdata1, ydata1, '%s' %(color1),label = 'Forcematch')   # create the initial plot
        xvals = np.arange(0,12,0.01)
	yvals = lj_a_coeff/xvals**12 - lj_b_coeff/xvals**6
	plt.plot(xvals,yvals,'%s' %(color2), label= 'LJ')

        plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')

	x_units = '$\AA$'
        y_units = 'kcal/mol'
        y_axis = 'Relative Free Energy'
        x_axis = 'Distance'
        x_axis = r'%s (%s)' %(x_axis, x_units)
        y_axis = r'%s (%s)' %(y_axis, y_units)
        plt.xlabel(r'%s' %(x_axis), size=12)
        plt.ylabel(r'%s' %(y_axis), size=12)

	plt.title(r'[%s] -- [%s]' %(name1,name2), size='14')

#       plt.xlim((0,200))
	plt.ylim((-10, 20.0))
	plt.legend(loc='upper right', ncol=1, fontsize = 'small')
       
        plt.savefig('lj_plots/%s_%s.lj.3.png' %(name1,name2))
        plt.close()



top_file = sys.argv[1]
cg_file = sys.argv[2]

ParsePrmtopBonded(top_file)

types = ['OXT','GCT','GCA','GHC','GN','GHN','GC','GO','PN','PCO','PO','PC1','PC2','PC3','PC4','PC5','PHC5','PHC4']
type_index_prmtop = [6,5,3,4,1,2,5,6,1,5,6,5,5,5,5,5,7,7]

if len(types) != len(type_index_prmtop):
	print "Types length (%s) != type_index_prmtop length (%s)" %(len(types),len(type_index_prmtop))

for i in range(len(types)):
	for j in range(i,len(types)):
		# Open cg_force_pair file and read in nonbonded interactions
		data = np.zeros((10000,2),dtype=float)
		print '[%s] -- [%s]' %(types[i],types[j])
		param2 = open(cg_file,'r')
		lines = param2.readlines()
		for k in range(len(lines)):
			n_letters = int(len(types[i])) + int(len(types[j])) + 13
			if lines[k][0:n_letters] == 'SF_Pairwise_%s_%s' %(types[i],types[j]):
				for m in range(10000):
					a,b,c,d = lines[k+3+m].split()
					data[m,0] = float(b)
					data[m,1] = float(c)
				
			elif lines[k][0:n_letters] == 'SF_Pairwise_%s_%s' %(types[j],types[i]):
				for m in range(10000):
					a,b,c,d = lines[k+3+m].split()
					data[m,0] = float(b)
					data[m,1] = float(c)

		if data[:,0].all() == 0:
			print 'error'
		# Get the index for the lennard jones parameters
		index = n_types*(type_index_prmtop[i]-1) + type_index_prmtop[j] - 1
		nb_index = nb_parm_index[index] - 1
		# Plot forcematched interactions and lj potenetials
		plot_1d(data[:,0],data[:,1],lj_a_coeff[nb_index],lj_b_coeff[nb_index],'k','r',types[i],types[j])
