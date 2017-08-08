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
from scipy import interpolate

kB = 0.001987 # kcal/mol/K
thresh = 1E-4

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, file_update, sim_in_file, atom_data_root, ib_lambda, n_iter, param_out_file, kT, psf_flag
	f = open(cfg_file)
	psf_flag = "true"
	file_update = ""
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			# check value
			if option.lower()=='topfile':
				top_file = value
			if option.lower()=='paramfile':
				param_out_file = value
			elif option.lower()=='siminfile':
				sim_in_file = value
			elif option.lower()=='atomisticdata':
				atom_data_root = value
			elif option.lower()=='lambda':
				ib_lambda = float(value)
			elif option.lower()=='iterations':
				n_iter = int(value)
			elif option.lower()=='psfflag':
				psf_flag = value
			elif option.lower()=='fileupdate':
				file_update = value
			elif option.lower()=='temperature':
				T = float(value)
			else :
				print "Option:", option, " is not recognized"
	f.close()
	if file_update != "true" and file_update != "false":
		file_update = "true"
		print "file_update being assigned default value of true"
	kT = kB*T

def ParseDataFiles():
	global uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds, uniq_angle_atom_types, n_uniq_angles, angles, n_angles,  uniq_dihedral_atom_types, n_uniq_dihedrals, dihedrals, n_dihedrals

	dihs = open("dihedrals_wuniq.txt",'r')
	line_count = 0
	for line in dihs:
		if line_count == 0 :	
			n_dihedrals = int(line[0:10])
			dihedrals = np.empty((n_dihedrals,5),dtype=int)
		else:
			for j in range(5):
				dihedrals[line_count-1,j] = int(line[j*11:j*11+10])
		line_count += 1
	dihs.close()
	dih_uniq = open("dihedrals_uniq_atom_types.txt",'r')
	line_count = 0
	uniq_dihedral_atom_types = []
	for line in dih_uniq:
		if line_count == 0 :	
			n_uniq_dihedrals = int(line[0:10])
		else:
			uniq_dihedral_atom_types.append([])
			for j in range(4):
				uniq_dihedral_atom_types[line_count-1].append(line[j*5:j*5+4].strip())
		line_count += 1
	dih_uniq.close()	

	angs = open("angles_wuniq.txt",'r')
	line_count = 0
	for line in angs:
		if line_count == 0 :	
			n_angles = int(line[0:10])
			angles = np.empty((n_angles,4),dtype=int)
		else:
			for j in range(4):
				angles[line_count-1,j] = int(line[j*11:j*11+10])
		line_count += 1
	angs.close()
	ang_uniq = open("angles_uniq_atom_types.txt",'r')
	line_count = 0
	uniq_angle_atom_types = []
	for line in ang_uniq:
		if line_count == 0 :	
			n_uniq_angles = int(line[0:10])
		else:
			uniq_angle_atom_types.append([])
			for j in range(3):
				uniq_angle_atom_types[line_count-1].append(line[j*5:j*5+4].strip())
		line_count += 1
	ang_uniq.close()

	bonds_file = open("bonds_wuniq.txt",'r')
	line_count = 0
	for line in bonds_file:
		if line_count == 0 :	
			n_bonds = int(line[0:10])
			bonds = np.empty((n_bonds,3),dtype=int)
		else:
			for j in range(3):
				bonds[line_count-1,j] = int(line[j*11:j*11+10])
		line_count += 1
	bonds_file.close()
	bond_uniq = open("bonds_uniq_atom_types.txt",'r')
	line_count = 0
	uniq_bond_atom_types = []
	for line in bond_uniq:
		if line_count == 0 :	
			n_uniq_bonds = int(line[0:10])
		else:
			uniq_bond_atom_types.append([])
			for j in range(2):
				uniq_bond_atom_types[line_count-1].append(line[j*5:j*5+4].strip())
		line_count += 1
	bond_uniq.close()


# read psf file and get bond, angle, and bond information
def ParsePsfFile(psf_file):
	global n_uniq_atom_types, atom_types, uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds, uniq_angle_atom_types, n_uniq_angles, angles, n_angles,  uniq_dihedral_atom_types, n_uniq_dihedrals, dihedrals, n_dihedrals
	f = open(psf_file)
	atom_flag = bond_flag = angle_flag = dihedral_flag = "false"
	atom = bond = angle = dihedral = 0
	bond_count = 0
	atom_types = []
	bonds = []
	uniq_bond_atom_types = []
	angles = []
	uniq_angle_atom_types = []
	dihedrals = []
	uniq_dihedral_atom_types = []
	n_uniq_atom_types = 0
	n_uniq_bonds = 0
	n_uniq_angles = 0
	n_uniq_dihedrals = 0
	for line in f:
		if atom_flag == "true" and atom < n_atoms:
			atom_types.append(line[29:33].strip())
			same = "false"
			if atom > 0:
				for atom2 in range(atom):
					if atom_types[atom] == atom_types[atom2]:
						same = "true"
						break
			if same == "false":
				n_uniq_atom_types += 1
			atom += 1
		elif bond_flag == "true" and bond < n_bonds:
			for line_bond in range(4):
				bonds.append([])
				bonds[bond].append(line[line_bond*16:line_bond*16+8])
				bonds[bond].append(line[line_bond*16+8:line_bond*16+16])
				same = "false"
				if bond > 0:
					for bond2 in range(bond):
						if (atom_types[int(bonds[bond][0])-1] == atom_types[int(bonds[bond2][0])-1] and atom_types[int(bonds[bond][1])-1] == atom_types[int(bonds[bond2][1])-1]) or (atom_types[int(bonds[bond][0])-1] == atom_types[int(bonds[bond2][1])-1] and atom_types[int(bonds[bond][1])-1] == atom_types[int(bonds[bond2][0])-1]):
							same = "true"
							uniq_bond_num = bonds[bond2][2]
							break
				if same == "false":
					uniq_bond_atom_types.append([])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][0])-1])
					uniq_bond_atom_types[n_uniq_bonds].append(atom_types[int(bonds[bond][1])-1])
					uniq_bond_num = str(n_uniq_bonds)
					n_uniq_bonds += 1
				bonds[bond].append(uniq_bond_num)
				bond += 1
				bond_count += 1
				if bond == n_bonds:
					break
		elif angle_flag == "true" and angle < n_angles:
			for line_angle in range(3):
				angles.append([])
				angles[angle].append(line[line_angle*24:line_angle*24+8])
				angles[angle].append(line[line_angle*24+8:line_angle*24+16])
				angles[angle].append(line[line_angle*24+16:line_angle*24+24])
				same = "false"
				if angle > 0:
					for angle2 in range(angle):
						if (atom_types[int(angles[angle][0])-1] == atom_types[int(angles[angle2][0])-1] and atom_types[int(angles[angle][1])-1] == atom_types[int(angles[angle2][1])-1] and atom_types[int(angles[angle][2])-1] == atom_types[int(angles[angle2][2])-1]) or (atom_types[int(angles[angle][0])-1] == atom_types[int(angles[angle2][2])-1] and atom_types[int(angles[angle][1])-1] == atom_types[int(angles[angle2][1])-1] and atom_types[int(angles[angle][2])-1] == atom_types[int(angles[angle2][0])-1]):
							same = "true"
							uniq_angle_num = angles[angle2][3]
							break
				if same == "false":
					uniq_angle_atom_types.append([])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][0])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][1])-1])
					uniq_angle_atom_types[n_uniq_angles].append(atom_types[int(angles[angle][2])-1])
					uniq_angle_num = str(n_uniq_angles)
					n_uniq_angles += 1
				angles[angle].append(uniq_angle_num)
				angle += 1
				if angle == n_angles:
					break
		elif dihedral_flag == "true" and dihedral < n_dihedrals:
			for line_dihedral in range(2):
				dihedrals.append([])
				dihedrals[dihedral].append(line[line_dihedral*32:line_dihedral*32+8])
				dihedrals[dihedral].append(line[line_dihedral*32+8:line_dihedral*32+16])
				dihedrals[dihedral].append(line[line_dihedral*32+16:line_dihedral*32+24])
				dihedrals[dihedral].append(line[line_dihedral*32+24:line_dihedral*32+32])
				same = "false"
				if dihedral > 0:
					for dihedral2 in range(dihedral):
						if (atom_types[int(dihedrals[dihedral][0])-1] == atom_types[int(dihedrals[dihedral2][0])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][3])-1]) or (atom_types[int(dihedrals[dihedral][0])-1] == atom_types[int(dihedrals[dihedral2][3])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][0])-1]):
							same = "true"
							uniq_dihedral_num = dihedrals[dihedral2][4]
							break
				if same == "false":
					uniq_dihedral_atom_types.append([])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][0])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][1])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][2])-1])
					uniq_dihedral_atom_types[n_uniq_dihedrals].append(atom_types[int(dihedrals[dihedral][3])-1])
					uniq_dihedral_num = str(n_uniq_dihedrals)
					n_uniq_dihedrals += 1
				dihedrals[dihedral].append(uniq_dihedral_num)
				dihedral += 1
				if dihedral == n_dihedrals:
					break

		if line[9:15] == "!NATOM":
			n_atoms = int(line[0:9])
			atom_flag = "true"
		elif line[9:15] == "!NBOND":
			n_bonds = int(line[0:9])
			bond_flag = "true"
		elif line[9:16] == "!NTHETA":
			n_angles = int(line[0:9])
			angle_flag = "true"
		elif line[9:14] == "!NPHI":
			n_dihedrals = int(line[0:9])
			dihedral_flag = "true"

	f.close()
	bonds = np.asmatrix(bonds,dtype=int)
	angles = np.asmatrix(angles,dtype=int)
	dihedrals = np.asmatrix(dihedrals,dtype=int)

	dih_temp = open("dihedrals_wuniq.txt",'w')
	dih_temp.write("%10d\n" %(n_dihedrals))
	for i in range(n_dihedrals):
		dih_temp.write("%10d %10d %10d %10d %10d\n" %(dihedrals[i,0], dihedrals[i,1], dihedrals[i,2], dihedrals[i,3], dihedrals[i,4]))
	dih_temp.close()
	dih_temp = open("dihedrals_uniq_atom_types.txt",'w')
	dih_temp.write("%10d\n" %(n_uniq_dihedrals))
	for i in range(n_uniq_dihedrals):
		dih_temp.write("%4s %4s %4s %4s\n" %(uniq_dihedral_atom_types[i][0], uniq_dihedral_atom_types[i][1], uniq_dihedral_atom_types[i][2], uniq_dihedral_atom_types[i][3]))
	dih_temp.close()
	ang_temp = open("angles_wuniq.txt",'w')
	ang_temp.write("%10d\n" %(n_angles))
	for i in range(n_angles):
		ang_temp.write("%10d %10d %10d %10d\n" %(angles[i,0], angles[i,1], angles[i,2], angles[i,3]))
	ang_temp.close()
	ang_temp = open("angles_uniq_atom_types.txt",'w')
	ang_temp.write("%10d\n" %(n_uniq_angles))
	for i in range(n_uniq_angles):
		ang_temp.write("%4s %4s %4s\n" %(uniq_angle_atom_types[i][0], uniq_angle_atom_types[i][1], uniq_angle_atom_types[i][2]))
	ang_temp.close()
	bond_temp = open("bonds_wuniq.txt",'w')
	bond_temp.write("%10d\n" % (n_bonds))
	for i in range(n_bonds):
		bond_temp.write("%10d %10d %10d\n" %(bonds[i,0], bonds[i,1], bonds[i,2]))
	bond_temp.close()
	bond_temp = open("bonds_uniq_atom_types.txt",'w')
	bond_temp.write("%10d\n" %(n_uniq_bonds))
	for i in range(n_uniq_bonds):
		bond_temp.write("%4s %4s\n" %(uniq_bond_atom_types[i][0], uniq_bond_atom_types[i][1]))
	bond_temp.close()

# add to dihedral distance histograms 
def ComputeDihedralHists(atom_positions, dihedrals, dihedral_hists, dihedral_min, delta_dihedral):
	# get size info
	n_dihedrals = dihedrals.shape[0]
	n_bins = dihedral_hists.shape[1]
	for dihedral in range(n_dihedrals):

		dih = ComputeDih(atom_positions[dihedrals[dihedral,0]-1,:],atom_positions[dihedrals[dihedral,1]-1,:],atom_positions[dihedrals[dihedral,2]-1,:],atom_positions[dihedrals[dihedral,3]-1,:])
		dih_bin = int((dih-dihedral_min)/delta_dihedral)
		if dih_bin >= 0 and dih_bin < n_bins:
			dihedral_hists[dihedrals[dihedral,4],dih_bin] += 1

# average dihedral distance histograms (convert them to probability densities)
def ReadDihedralHists(dihedral_hists, dihedral_min, delta_dihedral, uniq_dihedral_atom_types):
	global file_update

	n_dihedrals = dihedral_hists.shape[0]
	n_bins = dihedral_hists.shape[1]


	for dihedral in range(n_dihedrals):
		filename = "../../dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".hist"
		filename2 = "../../dihs/"+uniq_dihedral_atom_types[dihedral][3].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][0].strip()+".hist"
		if os.path.isfile(filename):
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_hists[dihedral,count] += float(val)
					count += 1
				inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					dihedral_hists[dihedral,count] += float(val)
					count += 1
				inp.close()
		else:
			print "Did not find data for dihedral", uniq_dihedral_atom_types[dihedral]
			sys.exit()


# add to angle distance histograms 
def ComputeAngleHists(atom_positions, angles, angle_hists, angle_min, delta_angle):
	# get size info
	n_angles = angles.shape[0]
	n_bins = angle_hists.shape[1]
	for angle in range(n_angles):

		ang = ComputeAng(atom_positions[angles[angle,0]-1,:],atom_positions[angles[angle,1]-1,:],atom_positions[angles[angle,2]-1,:])
		ang_bin = int((ang-angle_min)/delta_angle)
		if ang_bin >= 0 and ang_bin < n_bins:
			angle_hists[angles[angle,3],ang_bin] += 1

# average angle distance histograms (convert them to probability densities)
def ReadAngleHists(angle_hists, angle_min, delta_angle, uniq_angle_atom_types):
	global file_update

	n_angles = angle_hists.shape[0]
	n_bins = angle_hists.shape[1]


	for angle in range(n_angles):
		filename = "../../angs/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".hist"
		filename2 = "../../angs/"+uniq_angle_atom_types[angle][2].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][0].strip()+".hist"
		if os.path.isfile(filename):
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					angle_hists[angle,count] = float(val)
					count += 1
				inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					angle_hists[angle,count] = float(val)
					count += 1
				inp.close()
		else :
			print "Did not find data for angle", uniq_angle_atom_types[angle]
			sys.exit()

# add to bond distance histograms 
def ComputeBondHists(atom_positions, bonds, bond_hists, bond_min, delta_bond):
	# get sizes etc
	n_bonds = bonds.shape[0]
	n_bins = bond_hists.shape[1]
	for bond in range(n_bonds):

		dist = math.sqrt(ComputeDist2(atom_positions[bonds[bond,0]-1,:],atom_positions[bonds[bond,1]-1,:]))
		dist_bin = int((dist-bond_min)/delta_bond)
		if dist_bin >= 0 and dist_bin < n_bins:
			bond_hists[bonds[bond,2],dist_bin] += 1

# average bond distance histograms (convert them to probability densities)
def ReadBondHists(bond_hists, bond_min, delta_bond, uniq_bond_atom_types):
	global file_update

	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]

	for bond in range(n_bonds):
		filename = "../../bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".hist"
		filename2 = "../../bonds/"+uniq_bond_atom_types[bond][1].strip()+"_"+uniq_bond_atom_types[bond][0].strip()+".hist"
		if os.path.isfile(filename):
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] += float(val)
					count += 1
				inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			if file_update == "true":
				inp = open(filename,'r')
				count = 0
				for line in inp:
					junk, val = line.split()
					bond_hists[bond,count] += float(val)
					count += 1
				inp.close()
		else :
			print "Did not find prameters for:", uniq_bond_atom_types[bond]
			sys.exit()


# compute the dihedral between four points
def ComputeDih(r1,r2,r3,r4):

	r12 = np.empty(3,dtype=float)
	r23 = np.empty(3,dtype=float)
	r34 = np.empty(3,dtype=float)

	# compute two vectors (r1-r2 and r3-r2) making sure they are in same box
	for j in range(0,3):
		temp1 = r1[j]-r2[j]
		temp2 = r2[j]-r3[j]
		temp3 = r3[j]-r4[j]
		r12[j] = temp1
		r23[j] = temp2
		r34[j] = temp3
	
	A = np.cross(r12,r23)
	A /= math.sqrt(np.dot(A,A))
	B = np.cross(r23,r34)
	B /= math.sqrt(np.dot(B,B))
	C = np.cross(r23,A)
	C /= math.sqrt(np.dot(C,C))
	cos_phi = np.dot(A,B)
	sin_phi = np.dot(C,B)
	
	phi = -math.atan2(sin_phi,cos_phi)*180.0/3.1415926535
	return phi

# compute the angle between three points
def ComputeAng(r1,r2,r3):

	r21 = np.empty(3,dtype=float)
	r23 = np.empty(3,dtype=float)

	# compute two vectors (r1-r2 and r3-r2) making sure they are in same box
	for j in range(0,3):
		temp1 = r1[j]-r2[j]
		temp2 = r3[j]-r2[j]
		r21[j] = temp1
		r23[j] = temp2
	
	theta = math.acos(np.dot(r21,r23)/(math.sqrt(np.dot(r21,r21))*math.sqrt(np.dot(r23,r23))))*180.0/3.1415926535
	return theta

# compute the distance between two points taking into account periodic boundary conditions
def ComputeDist2(r1,r2):
	dist2 = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		dist2 += temp*temp

	return dist2;

# compute new CG probability distributions from CG trajectory
def AnalyzeCGTraj(top_file, traj_file, cg_bond_hists, cg_angle_hists, cg_dihedral_hists):
	global bonds, bond_min, delta_bond, angles, angle_min, delta_angle, dihedrals, dihedral_min, delta_dihedral
	# initiate MDAnalysis coordinate universe
	coord = MDAnalysis.Universe(top_file, traj_file)
	# make an atom selection
	sel = coord.select_atoms("all")
	n_frames = coord.trajectory.n_frames
	print "Number of frames in trajectory file: ", n_frames
	# declare array that will contain the selected atom positions at each time frame
	positions = np.empty((sel.n_atoms,3),dtype=float)
	# zero hist array
	cg_bond_hists.fill(0.0)
	cg_angle_hists.fill(0.0)
	cg_dihedral_hists.fill(0.0)

	# Loop through trajectory
	for ts in coord.trajectory:
	
		# 
		if ts.frame % 1000 == 0:
			print "Frame: ", ts.frame, " of ", n_frames
		
		# save selected positions into array
		positions = sel.positions
		
		# add to bond histograms
		ComputeBondHists(positions, bonds, cg_bond_hists, bond_min, delta_bond)
		# add to angle histograms
		ComputeAngleHists(positions, angles, cg_angle_hists, angle_min, delta_angle)
		# add to dihedral histograms
		ComputeDihedralHists(positions, dihedrals, cg_dihedral_hists, dihedral_min, delta_dihedral)

def ReadAtomHists(atom_bond_hists, bond_potentials, atom_angle_hists, angle_potentials, atom_dihedral_hists, dihedral_potentials):
	global bond_min, delta_bond, uniq_bond_atom_types, angle_min, delta_angle, uniq_angle_atom_types, dihedral_min, delta_dihedral, uniq_dihedral_atom_types, param_out_file

	ReadBondHists(atom_bond_hists, bond_min, delta_bond, uniq_bond_atom_types)
	ReadAngleHists(atom_angle_hists, angle_min, delta_angle, uniq_angle_atom_types)
	ReadDihedralHists(atom_dihedral_hists, dihedral_min, delta_dihedral, uniq_dihedral_atom_types)

	param_out = open("GGGGG.30.cg.bond1.lammpspar", 'w')
	CreateBondPotentials(atom_bond_hists, bond_potentials, bond_min, delta_bond, uniq_bond_atom_types,"bond0.ib","bond0.ener",param_out)
	CreateAnglePotentials(atom_angle_hists, angle_potentials, angle_min, delta_angle, uniq_angle_atom_types,"angle0.ib","angle0.ener",param_out)
	CreateDihedralPotentials(atom_dihedral_hists, dihedral_potentials, dihedral_min, delta_dihedral, uniq_dihedral_atom_types,"dihedral0.ib","dihedral0.ener",param_out)
	param_out.close()


# average dihedral distance histograms (convert them to probability densities)
def CreateDihedralPotentials(dihedral_hists, dihedral_potentials, dihedral_min, delta_dihedral, uniq_dihedral_atom_types, ib_out, ener_out, param_out):
	global file_update, kT

	n_dihedrals = dihedral_hists.shape[0]
	n_bins = dihedral_hists.shape[1]
	ener_temp = np.empty((n_dihedrals,n_bins),dtype=float)

	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty(n_bins,dtype=float)

	out = open(ib_out,'w')
	# first check to see if we have two copies of same dihedral
	for dihedral in range(n_dihedrals):
		for i in range(n_bins):
			x = dihedral_min+(i+0.5)*delta_dihedral
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			x_mat[i] = x
			coeff_mat[i,2] = x*x
		# write title to output 
#		out.write("\n%5d\n" % (dihedral+1))
		out.write("\nDIH_%s\n" % (str(dihedral+1)))
		out.write("N %5d DEGREES\n\n" % (n_bins))

		# convert to probability density
		dihedral_hists[dihedral,:] /= (np.sum(dihedral_hists[dihedral,:])*delta_dihedral)
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if dihedral_hists[dihedral,i] > thresh:
#				dihedral_hists[dihedral,i] = -kT*math.log(dihedral_hists[dihedral,i]/math.sin(coeff_mat[i,1]*3.1415926535/180.0))
				ener_temp[dihedral,i] = -kT*math.log(dihedral_hists[dihedral,i])
				if count == 0 or ener_temp[dihedral,i] < min_val:
					min_val = ener_temp[dihedral,i]
				count += 1
			else:
				ener_temp[dihedral,i] = 0
			# make sure the tails are increasing
			if i > 0 and ener_temp[dihedral,i] == 0 and ener_temp[dihedral,i-1] != 0:
				for j in range(1,4):
					if ener_temp[dihedral,i-j] < ener_temp[dihedral,i-1-j]:
						ener_temp[dihedral,i-j] = 0
#		ener_temp[dihedral,:] -= min_val
#		dihedral_potentials[dihedral,:] -= min_val
		j = 0
		while ener_temp[dihedral,j] == 0:
			j += 1
		first_nonzero = j
		for j in range(first_nonzero,first_nonzero+3):
			if ener_temp[dihedral,j+1] > ener_temp[dihedral,j]:
				ener_temp[dihedral,j] = 0
				dihedral_potentials[dihedral,j] = 0


		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[dihedral,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[dihedral,i] != 0:
				stop = i+1
			
			# fit previous section if we have enough data
			if ener_temp[dihedral,i] == 0 and data_flag == "true" and (stop-start) > 5:
				# smooth region with cubic spline interpolation
				smooth = 2
				tck = interpolate.splrep(x_mat[start:stop],ener_temp[dihedral,start:stop],s=smooth)
				ener_temp[dihedral,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
				# add data to missing values:
				if data_count == 0:
					init_points = ener_temp[dihedral,start:start+5]
					init_x = coeff_mat[start:start+5]
					init_start = start
					previous_points = ener_temp[dihedral,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop
				else:
					# fit to harmonic if big enough gap
					if (start-previous_stop) > 6:
						current_points = ener_temp[dihedral,start:start+5]
						current_x = coeff_mat[start:start+5]
						data_y = np.append(previous_points,current_points)
						data_x = np.append(previous_x, current_x, axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
#							force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[dihedral,previous_stop-2:start+2],s=smooth)
						ener_temp[dihedral,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
						force[previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=1)
					# otherwise fit to more steep quadratic?
					else:
						current_points = ener_temp[dihedral,start:start+2]
						current_x = coeff_mat[start:start+2]
						data_y = np.append(previous_points[-2:],current_points)
						data_x = np.append(previous_x[-2:],current_x,axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
#							force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[dihedral,previous_stop-2:start+2],s=smooth)
						ener_temp[dihedral,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
						force[previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=1)
					# save current data into previous
					previous_points = ener_temp[dihedral,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop

				# compute force using finite difference
				for j in range(start,stop):
					if j > start and j < stop-1:
						force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j-1])/(2*delta_dihedral)
					elif j < stop-1:
						force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j])/(delta_dihedral)
					else:
						force[j] = -(ener_temp[dihedral,j]-ener_temp[dihedral,j-1])/(delta_dihedral)

				data_flag = "false"
				data_count += 1
			# skip 'island' of data
			elif ener_temp[dihedral,i] == 0 and data_flag == "true":
				data_flag = "false"

			i += 1

		if data_flag == "true" and (n_bins-start) > 5:
			init_points = ener_temp[dihedral,start:start+5]
			init_x = coeff_mat[start:start+5]
			init_start = start
			previous_stop = n_bins
			previous_x = coeff_mat[n_bins-5:n_bins]
			previous_points = ener_temp[dihedral,n_bins-5:n_bins]

		# fit remaining missing data
		for j in range(5):
			previous_x[j,1] -= 360.0
			previous_x[j,2] = previous_x[j,1]*previous_x[j,1]
		data_y = np.append(previous_points,init_points)
		data_x = np.append(previous_x, init_x, axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(init_start):
			ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
			force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
		for j in range(previous_stop,n_bins):
			x = dihedral_min+(j+0.5)*delta_dihedral-360.0
			x2 = x*x
			ener_temp[dihedral,j] = k[0] + k[1]*x + k[2]*x2
			force[j] = -2*k[2]*x - k[1]

		# now smooth the potential using a cubic spline fit	
		smooth=3
		# first connect point
		if init_start > 4:
			tck = interpolate.splrep(x_mat[init_start-4:init_start+4],ener_temp[dihedral,init_start-4:init_start+4],s=smooth)
			ener_temp[dihedral,init_start-4:init_start+4] = interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=0)
			force[init_start-4:init_start+4] = -interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=1)
		else:
			tck = interpolate.splrep(x_mat[0:init_start+4],ener_temp[dihedral,0:init_start+4],s=smooth)
			ener_temp[dihedral,0:init_start+4] = interpolate.splev(x_mat[0:init_start+4], tck, der=0)
			force[0:init_start+4] = -interpolate.splev(x_mat[0:init_start+4], tck, der=1)
		# last connect point
		if previous_stop < n_bins - 4 :
			tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[dihedral,previous_stop-4:previous_stop+4],s=smooth)
			ener_temp[dihedral,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
			force[previous_stop-4:previous_stop+4] = -interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=1)
		else:
			tck = interpolate.splrep(x_mat[previous_stop-4:n_bins],ener_temp[dihedral,previous_stop-4:n_bins],s=smooth)
			ener_temp[dihedral,previous_stop-4:n_bins] = interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=0)
			force[previous_stop-4:n_bins] = -interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=1)
		
		min_val = np.min(ener_temp[dihedral,:])
		ener_temp[dihedral,:] -= min_val
                # compute force using finite difference
		for j in range(n_bins):
                        if j > 0 and j < n_bins-1:
                                force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j-1])/(2*delta_dihedral)
                        elif j == 0:
                                force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j])/(delta_dihedral)
                        else:
                                force[j] = -(ener_temp[dihedral,j]-ener_temp[dihedral,j-1])/(delta_dihedral)

		for i in range(n_bins):
			dihedral_potentials[dihedral,i] = ener_temp[dihedral, i]
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,dihedral_min+(i)*delta_dihedral, ener_temp[dihedral,i],force[i]))
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, ib_out, str(dihedral+1).strip()))

		
	out.close()
	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (dihedral_min+(i+0.5)*delta_dihedral))
		for dihedral in range(n_dihedrals):
			temp.write("%12.7f" % (ener_temp[dihedral,i]))
		temp.write("\n")
	temp.close()
# MM_DIHEDRALS

# average angle distance histograms (convert them to probability densities)
def CreateAnglePotentials(angle_hists, angle_potentials, angle_min, delta_angle, uniq_angle_atom_types, ib_out, ener_out, param_out):
	global file_update, kT

	n_angles = angle_hists.shape[0]
	n_bins = angle_hists.shape[1]
	ener_temp = np.empty((n_angles,n_bins),dtype=float)

	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty((n_bins),dtype=float)
	out = open(ib_out,'w')
	# first check to see if we have two copies of same angle
	for angle in range(n_angles):
		for i in range(n_bins):
			x = angle_min+(i+0.5)*delta_angle
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			x_mat[i] = x
			coeff_mat[i,2] = x*x
		out.write("\nANG_%s\n" % (str(angle+1)))
		out.write("N %5d\n\n" % (n_bins))

		# convert to probability density
		angle_hists[angle,:] /= (np.sum(angle_hists[angle,:])*delta_angle)
#		ener_temp[angle,:] = angle_hists[angle,:]
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if angle_hists[angle,i] > thresh:
#				angle_hists[angle,i] = -kT*math.log(angle_hists[angle,i]/math.sin(coeff_mat[i,1]*3.1415926535/180.0))
				ener_temp[angle,i] = -kT*math.log(angle_hists[angle,i]/math.sin(coeff_mat[i,1]*np.pi/180))
				angle_potentials[angle,i] = ener_temp[angle,i]
				if count == 0 or ener_temp[angle,i] < min_val:
					min_val = ener_temp[angle,i]
				count += 1
			else:
				ener_temp[angle,i] = 0
				angle_potentials[angle,i] = 0
			# make sure the tails are increasing
			if i > 0 and ener_temp[angle,i] == 0 and ener_temp[angle,i-1] != 0:
				for j in range(1,4):
					if ener_temp[angle,i-j] < ener_temp[angle,i-1-j]:
						ener_temp[angle,i-j] = 0
						angle_potentials[angle,i-j] = 0
		# zero potentials
#		ener_temp[angle,:] -= min_val
#		angle_potentials[angle,:] -= min_val
	
		j = 0
		while ener_temp[angle,j] == 0:
			j += 1
		first_nonzero = j
		for j in range(first_nonzero,first_nonzero+3):
			if ener_temp[angle,j+1] > ener_temp[angle,j]:
				ener_temp[angle,j] = 0
				angle_potentials[angle,j] = 0


		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[angle,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[angle,i] != 0:
				stop = i+1
			
			# fit previous section if we have enough data
			if ener_temp[angle,i] == 0 and data_flag == "true" and (stop-start) > 5:
				# smooth region with cubic spline interpolation
				smooth = 2
				tck = interpolate.splrep(x_mat[start:stop],ener_temp[angle,start:stop],s=smooth)
				ener_temp[angle,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
#				k, rss, rank, s = np.linalg.lstsq(coeff_mat[start:stop],ener_temp[angle,start:stop])
				# add data to missing values:
				if data_count == 0:
					init_points = ener_temp[angle,start:start+5]
					init_x = coeff_mat[start:start+5]
					init_start = start
					previous_points = ener_temp[angle,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop
				else:
					# fit to harmonic if big enough gap
					if (start-previous_stop) > 6:
						current_points = ener_temp[angle,start:start+5]
						current_x = coeff_mat[start:start+5]
						data_y = np.append(previous_points,current_points)
						data_x = np.append(previous_x, current_x, axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[angle,previous_stop-2:start+2],s=smooth)
						ener_temp[angle,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
						force[previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=1)
#						force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
					# otherwise fit to more steep quadratic?
					else:
						current_points = ener_temp[angle,start:start+2]
						current_x = coeff_mat[start:start+2]
						data_y = np.append(previous_points[-2:],current_points)
						data_x = np.append(previous_x[-2:],current_x,axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
#							force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[angle,previous_stop-2:start+2],s=smooth)
						ener_temp[angle,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
						force[previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=1)
					# save current data into previous
					previous_points = ener_temp[angle,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop

				# compute force using finite difference
				for j in range(start,stop):
					if j > start and j < stop-1:
						force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j-1])/(2*delta_angle)
					elif j < stop-1:
						force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j])/(delta_angle)
					else:
						force[j] = -(ener_temp[angle,j]-ener_temp[angle,j-1])/(delta_angle)

				data_flag = "false"
				data_count += 1
			# skip 'island' of data
			elif ener_temp[angle,i] == 0 and data_flag == "true":
				data_flag = "false"

			i += 1

		# check to see if we are at 180 and thus need to ensure periodicity
		if data_flag == "true" and (n_bins-start) > 5:
			init_points = ener_temp[angle,start:start+5]
			init_x = coeff_mat[start:start+5]
			init_start = start
			previous_stop = n_bins
			previous_x = coeff_mat[n_bins-5:n_bins]
			previous_points = ener_temp[angle,n_bins-5:n_bins]

		# fit remaining missing data
		for j in range(5):
			previous_x[j,1] -= 180.0
			previous_x[j,2] = previous_x[j,1]*previous_x[j,1]
		data_y = np.append(previous_points,init_points)
		data_x = np.append(previous_x, init_x, axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(init_start):
			ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
			force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
		for j in range(previous_stop,n_bins):
			x = angle_min+(j+0.5)*delta_angle-180.0
			x2 = x*x
			ener_temp[angle,j] = k[0] + k[1]*x + k[2]*x2
			force[j] = -2*k[2]*x - k[1]

		# now smooth the potential using a cubic spline fit	
		smooth=3
		# first connect point
		if init_start > 4:
			tck = interpolate.splrep(x_mat[init_start-4:init_start+4],ener_temp[angle,init_start-4:init_start+4],s=smooth)
			ener_temp[angle,init_start-4:init_start+4] = interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=0)
			force[init_start-4:init_start+4] = -interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=1)
		else:
			tck = interpolate.splrep(x_mat[0:init_start+4],ener_temp[angle,0:init_start+4],s=smooth)
			ener_temp[angle,0:init_start+4] = interpolate.splev(x_mat[0:init_start+4], tck, der=0)
			force[0:init_start+4] = -interpolate.splev(x_mat[0:init_start+4], tck, der=1)
		# last connect point
		if previous_stop < n_bins - 4 :
			tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[angle,previous_stop-4:previous_stop+4],s=smooth)
			ener_temp[angle,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
			force[previous_stop-4:previous_stop+4] = -interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=1)
		else:
			tck = interpolate.splrep(x_mat[previous_stop-4:n_bins],ener_temp[angle,previous_stop-4:n_bins],s=smooth)
			ener_temp[angle,previous_stop-4:n_bins] = interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=0)
			force[previous_stop-4:n_bins] = -interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=1)

		min_val = np.min(ener_temp[angle,:])
		ener_temp[angle,:] -= min_val
		# compute force using finite difference
		for j in range(n_bins):
                        if j > 0 and j < n_bins-1:
                                force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j-1])/(2*delta_angle)
                        elif j == 0:
                                force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j])/(delta_angle)
                        else:
                                force[j] = -(ener_temp[angle,j]-ener_temp[angle,j-1])/(delta_angle)
		
		for i in range(n_bins):
			angle_potentials[angle,i] = ener_temp[angle,i]
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,angle_min+(i)*delta_angle, ener_temp[angle,i],force[i]))
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, ib_out, str(angle+1).strip()))
	out.close()
	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (angle_min+(i+0.5)*delta_angle))
		for angle in range(n_angles):
			temp.write("%12.7f" % (ener_temp[angle,i]))
		temp.write("\n")
	temp.close()

# MM_ANGLES

# average bond distance histograms (convert them to probability densities)
def CreateBondPotentials(bond_hists, bond_potentials, bond_min, delta_bond, uniq_bond_atom_types,ib_out,ener_out,param_out):
	global file_update, kT

	n_bonds = bond_hists.shape[0]
	n_bins = bond_hists.shape[1]
	ener_temp = np.empty((n_bonds,n_bins),dtype=float)
	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty(n_bond_bins,dtype=float)
	
	out = open(ib_out,'w')
	for bond in range(n_bonds):
		for i in range(n_bins):
			x = bond_min+(i+0.5)*delta_bond
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			coeff_mat[i,2] = x*x
			x_mat[i] = x
		# write title to output 
		out.write("\nBOND_%s\n" % (str(bond+1)))
		out.write("N %5d\n\n" % (n_bins))
		# create probability density
		bond_hists[bond,:] /= (np.sum(bond_hists[bond,:])*delta_bond)
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if bond_hists[bond,i] > thresh :
				ener_temp[bond,i] = -kT*math.log(bond_hists[bond,i]/coeff_mat[i,2])
				if count == 0 or ener_temp[bond,i] < min_val:
					min_val = ener_temp[bond,i]
				count += 1
			else: 
				ener_temp[bond,i] = 0
		# zero potential
#		ener_temp[bond,:] -= min_val
#		bond_potentials[bond,:] -= min_val
		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[bond,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[bond,i] != 0:
				stop = i+1
			
                	# fit previous section if we have enough data
			if ener_temp[bond,i] == 0 and data_flag == "true" and (stop-start) > 4:

				# add data to missing values
				k, rss, rank, s = np.linalg.lstsq(coeff_mat[start:stop],ener_temp[bond,start:stop])
				print k
				for j in range(start):
					ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
					force[j] = -2.0*k[2]*coeff_mat[j,1] - k[1]
				for j in range(stop,n_bins):
					ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
					force[j] = -2.0*k[2]*coeff_mat[j,1] - k[1]
			
				# compute force
				for j in range(start,stop):
					if j > start and j < stop-1:
						force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j-1])/(2*delta_bond)
					elif j == start:
						force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j])/(delta_bond)
					elif j == stop-1:
						force[j] = -(ener_temp[bond,j]-ener_temp[bond,j-1])/(delta_bond)
			
				data_flag = "false"
				data_count += 1
				break
			# skip 'island' of data
			elif ener_temp[bond,i] == 0 and data_flag == "true":
				data_flag = "false"
			
			i += 1
	
		# now smooth the potential using a cubic spline fit	
		tck = interpolate.splrep(coeff_mat[:,1],ener_temp[bond,:]-min_val,s=n_bins-math.sqrt(2*n_bins))
		f2 = interpolate.splev(coeff_mat[:,1], tck, der=0)
		force = -interpolate.splev(coeff_mat[:,1], tck, der=1)
		for i in range(n_bins):
			ener_temp[bond,i] = f2[i]
		min_val = np.min(ener_temp[bond,:])
		ener_temp[bond,:] -= min_val
		
		for i in range(n_bins):
			bond_potentials[bond,i] = ener_temp[bond,i] 
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*delta_bond, ener_temp[bond,i],force[i]))		
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, ib_out, str(bond+1).strip()))
	out.close()

	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (bond_min+(i+0.5)*delta_bond))
		for bond in range(n_bonds):
			temp.write(" %12.7f" % (ener_temp[bond,i]))
		temp.write("\n")
	temp.close()

# update dihedral potentials using iterative boltzmann procedure
def UpdateDihedrals(ib_iter, atom_dihedral_potentials, cg_dihedral_hists, cg_dihedral_potentials, param_out):
	global file_update, kT, ib_lambda, dihedral_min, delta_dihedral

	n_dihedrals = cg_dihedral_hists.shape[0]
	n_bins = cg_dihedral_hists.shape[1]
	ener_temp = np.zeros((n_dihedrals,n_bins),dtype=float)

	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "dihedral" + str(ib_iter+1) + ".ib"
	ener_out = "dihedral" + str(ib_iter+1) + ".ener"
	out = open(ib_out,'w')
	# first check to see if we have two copies of same dihedral
	for dihedral in range(n_dihedrals):
		# convert to probability density
		cg_dihedral_hists[dihedral,:] /= (np.sum(cg_dihedral_hists[dihedral,:])*delta_dihedral)
		for i in range(n_bins):
			x = dihedral_min+(i+0.5)*delta_dihedral
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			x_mat[i] = x
			coeff_mat[i,2] = x*x
		# write title to output 
		out.write("\nDIH_%s\n" % (str(dihedral+1)))
		out.write("N %5d DEGREES\n\n" % (n_bins))

		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if cg_dihedral_hists[dihedral,i] > thresh:
				ener_temp[dihedral,i] = -kT*math.log(cg_dihedral_hists[dihedral,i])
			else:
				ener_temp[dihedral,i] = 0
			# make sure the tails are increasing
			if i > 0 and ener_temp[dihedral,i] == 0 and ener_temp[dihedral,i-1] != 0:
				for j in range(1,4):
					if ener_temp[dihedral,i-j] < ener_temp[dihedral,i-1-j]:
						ener_temp[dihedral,i-j] = 0

		j = 0
		while ener_temp[dihedral,j] == 0:
			j += 1
		first_nonzero = j
		for j in range(first_nonzero,first_nonzero+3):
			if ener_temp[dihedral,j+1] > ener_temp[dihedral,j]:
				ener_temp[dihedral,j] = 0


		# zero potential
		min_val = np.amin(ener_temp[dihedral,:])

		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[dihedral,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[dihedral,i] != 0:
				stop = i+1
			
			# fit previous section if we have enough data
			if ener_temp[dihedral,i] == 0 and data_flag == "true" and (stop-start) > 5:
				# smooth region with cubic spline interpolation
				smooth = 10
				tck = interpolate.splrep(x_mat[start:stop],ener_temp[dihedral,start:stop],s=smooth)
				ener_temp[dihedral,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
				# add data to missing values:
				if data_count == 0:
					init_points = ener_temp[dihedral,start:start+5]
					init_x = coeff_mat[start:start+5]
					init_start = start
					previous_points = ener_temp[dihedral,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop
				else:
					# fit to harmonic if big enough gap
					if (start-previous_stop) > 6:
						current_points = ener_temp[dihedral,start:start+5]
						current_x = coeff_mat[start:start+5]
						data_y = np.append(previous_points,current_points)
						data_x = np.append(previous_x, current_x, axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-3:start+3],ener_temp[dihedral,previous_stop-3:start+3],s=smooth)
						ener_temp[dihedral,previous_stop-3:start+3] = interpolate.splev(x_mat[previous_stop-3:start+3], tck, der=0)
					# otherwise fit to more steep quadratic?
					else:
						current_points = ener_temp[dihedral,start:start+2]
						current_x = coeff_mat[start:start+2]
						data_y = np.append(previous_points[-2:],current_points)
						data_x = np.append(previous_x[-2:],current_x,axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[dihedral,previous_stop-2:start+2],s=smooth)
						ener_temp[dihedral,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
					# save current data into previous
					previous_points = ener_temp[dihedral,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop

				data_flag = "false"
				data_count += 1
			# skip 'island' of data
			elif ener_temp[dihedral,i] == 0 and data_flag == "true":
				data_flag = "false"

			i += 1

		if data_flag == "true" and (n_bins-start) > 5:
			init_points = ener_temp[dihedral,start:start+5]
			init_x = coeff_mat[start:start+5]
			init_start = start
			previous_stop = n_bins
			previous_x = coeff_mat[n_bins-5:n_bins]
			previous_points = ener_temp[dihedral,n_bins-5:n_bins]

		# fit remaining missing data
		for j in range(5):
			previous_x[j,1] -= 360.0
			previous_x[j,2] = previous_x[j,1]*previous_x[j,1]
		data_y = np.append(previous_points,init_points)
		data_x = np.append(previous_x, init_x, axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(init_start):
			ener_temp[dihedral,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		for j in range(previous_stop,n_bins):
			x = dihedral_min+(j+0.5)*delta_dihedral-360.0
			x2 = x*x
			ener_temp[dihedral,j] = k[0] + k[1]*x + k[2]*x2

		# now smooth the potential using a cubic spline fit	
		smooth=3
		# first connect point
		if init_start > 4:
			tck = interpolate.splrep(x_mat[init_start-4:init_start+4],ener_temp[dihedral,init_start-4:init_start+4],s=smooth)
			ener_temp[dihedral,init_start-4:init_start+4] = interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=0)
		else:
			tck = interpolate.splrep(x_mat[0:init_start+4],ener_temp[dihedral,0:init_start+4],s=smooth)
			ener_temp[dihedral,0:init_start+4] = interpolate.splev(x_mat[0:init_start+4], tck, der=0)
		# last connect point
		if previous_stop < n_bins - 4 :
			tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[dihedral,previous_stop-4:previous_stop+4],s=smooth)
			ener_temp[dihedral,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
		else:
			tck = interpolate.splrep(x_mat[previous_stop-4:n_bins],ener_temp[dihedral,previous_stop-4:n_bins],s=smooth)
			ener_temp[dihedral,previous_stop-4:n_bins] = interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=0)

		min_val = np.min(ener_temp[dihedral,:])
                ener_temp[dihedral,:] -= min_val
		for i in range(n_bins):
			cg_dihedral_potentials[dihedral,i] -= ib_lambda * (ener_temp[dihedral,i] - atom_dihedral_potentials[dihedral,i])
		min_val = np.min(cg_dihedral_potentials[dihedral,:])
                cg_dihedral_potentials[dihedral,:] -= min_val
		for i in range(n_bins):
			# finite difference force:
			if i==0:
				force[i] = -(cg_dihedral_potentials[dihedral,i+1]-cg_dihedral_potentials[dihedral,i])/delta_dihedral
			elif i == n_bins-1:
				force[i] = -(cg_dihedral_potentials[dihedral,i]-cg_dihedral_potentials[dihedral,i-1])/delta_dihedral
			else:
				force[i] = -(cg_dihedral_potentials[dihedral,i+1]-cg_dihedral_potentials[dihedral,i-1])/(2.0*delta_dihedral)
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,dihedral_min+(i)*delta_dihedral, cg_dihedral_potentials[dihedral,i],force[i]))
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, ib_out, str(dihedral+1).strip()))
		out.flush()
	out.close()
	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (dihedral_min+(i+0.5)*delta_dihedral))
		for dihedral in range(n_dihedrals):
			temp.write("%12.7f" % (cg_dihedral_potentials[dihedral,i]))
		temp.write("\n")
	temp.close()
# MM_DIHEDRALS

def UpdateAngles(ib_iter, atom_angle_potentials, cg_angle_hists, cg_angle_potentials, param_out):
	global file_update, kT, ib_lambda, angle_min, delta_angle

	n_angles = cg_angle_hists.shape[0]
	n_bins = cg_angle_hists.shape[1]
	ener_temp = np.zeros((n_angles,n_bins),dtype=float)

	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "angle" + str(ib_iter+1) + ".ib" 
	ener_out = "angle" + str(ib_iter+1) + ".ener" 
	out = open(ib_out,'w')
	# first check to see if we have two copies of same angle
	for angle in range(n_angles):
		# convert to probability density
		cg_angle_hists[angle,:] /= (np.sum(cg_angle_hists[angle,:])*delta_angle)
		for i in range(n_bins):
			x = angle_min+(i+0.5)*delta_angle
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			x_mat[i] = x
			coeff_mat[i,2] = x*x
		out.write("\nANG_%s\n" % (str(angle+1)))
		out.write("N %5d\n\n" % (n_bins))

		# convert to energies
		for i in range(n_bins):
			if cg_angle_hists[angle,i] > thresh:
				ener_temp[angle,i] = -kT*math.log(cg_angle_hists[angle,i])
			else:
				ener_temp[angle,i] = 0
			# make sure the tails are increasing
			if i > 0 and ener_temp[angle,i] == 0 and ener_temp[angle,i-1] != 0:
				for j in range(1,4):
					if ener_temp[angle,i-j] < ener_temp[angle,i-1-j]:
						ener_temp[angle,i-j] = 0

		# zero the minimum
		min_val = np.amin(ener_temp[angle,:])
		j = 0
		while ener_temp[angle,j] == 0:
			j += 1
		first_nonzero = j
		for j in range(first_nonzero,first_nonzero+3):
			if ener_temp[angle,j+1] > ener_temp[angle,j]:
				ener_temp[angle,j] = 0

		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[angle,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[angle,i] != 0:
				stop = i+1
			
			# fit previous section if we have enough data
			if ener_temp[angle,i] == 0 and data_flag == "true" and (stop-start) > 5:
				# smooth region with cubic spline interpolation
				smooth = 10
				tck = interpolate.splrep(x_mat[start:stop],ener_temp[angle,start:stop],s=smooth)
				ener_temp[angle,start:stop] = interpolate.splev(x_mat[start:stop], tck, der=0)
				# add data to missing values:
				if data_count == 0:
					init_points = ener_temp[angle,start:start+5]
					init_x = coeff_mat[start:start+5]
					init_start = start
					previous_points = ener_temp[angle,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop
				else:
					# fit to harmonic if big enough gap
					if (start-previous_stop) > 6:
						current_points = ener_temp[angle,start:start+5]
						current_x = coeff_mat[start:start+5]
						data_y = np.append(previous_points,current_points)
						data_x = np.append(previous_x, current_x, axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-3:start+3],ener_temp[angle,previous_stop-3:start+3],s=smooth)
						ener_temp[angle,previous_stop-3:start+3] = interpolate.splev(x_mat[previous_stop-3:start+3], tck, der=0)
					# otherwise fit to more steep quadratic?
					else:
						current_points = ener_temp[angle,start:start+2]
						current_x = coeff_mat[start:start+2]
						data_y = np.append(previous_points[-2:],current_points)
						data_x = np.append(previous_x[-2:],current_x,axis=0)
						k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
						# fill in hole
						for j in range(previous_stop,start):
							ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
						# smooth region with cubic spline interpolation
						smooth = 2
						tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[angle,previous_stop-2:start+2],s=smooth)
						ener_temp[angle,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
					# save current data into previous
					previous_points = ener_temp[angle,stop-5:stop]
					previous_x = coeff_mat[stop-5:stop]
					previous_stop = stop

				data_flag = "false"
				data_count += 1
			# skip 'island' of data
			elif ener_temp[angle,i] == 0 and data_flag == "true":
				data_flag = "false"

			i += 1

		if data_flag == "true" and (n_bins-start) > 5:
			init_points = ener_temp[angle,start:start+5]
			init_x = coeff_mat[start:start+5]
			init_start = start
			previous_stop = n_bins
			previous_x = coeff_mat[n_bins-5:n_bins]
			previous_points = ener_temp[angle,n_bins-5:n_bins]
		elif data_flag=="true":
			print "oops ", angle
			print "start:", start
			print "n_bins:", n_bins
		# fit remaining missing data
		for j in range(5):
			previous_x[j,1] -= 180.0
			previous_x[j,2] = previous_x[j,1]*previous_x[j,1]
		data_y = np.append(previous_points,init_points)
		data_x = np.append(previous_x, init_x, axis=0)
		k, rss, rank, s = np.linalg.lstsq(data_x, data_y)

		for j in range(init_start):
			ener_temp[angle,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
		for j in range(previous_stop,n_bins):
			x = angle_min+(j+0.5)*delta_angle-180.0
			x2 = x*x
			ener_temp[angle,j] = k[0] + k[1]*x + k[2]*x2
		# now smooth the potential using a cubic spline fit	
		smooth=3
		# first connect point
		if init_start > 4:
			tck = interpolate.splrep(x_mat[init_start-4:init_start+4],ener_temp[angle,init_start-4:init_start+4],s=smooth)
			ener_temp[angle,init_start-4:init_start+4] = interpolate.splev(x_mat[init_start-4:init_start+4], tck, der=0)
		else:
			tck = interpolate.splrep(x_mat[0:init_start+4],ener_temp[angle,0:init_start+4],s=smooth)
			ener_temp[angle,0:init_start+4] = interpolate.splev(x_mat[0:init_start+4], tck, der=0)
		# last connect point
		if previous_stop < n_bins - 4 :
			tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[angle,previous_stop-4:previous_stop+4],s=smooth)
			ener_temp[angle,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
		else:
			tck = interpolate.splrep(x_mat[previous_stop-4:n_bins],ener_temp[angle,previous_stop-4:n_bins],s=smooth)
			ener_temp[angle,previous_stop-4:n_bins] = interpolate.splev(x_mat[previous_stop-4:n_bins], tck, der=0)
		
		min_val = np.min(ener_temp[angle,:])
                ener_temp[angle,:] -= min_val

		# write output
		for i in range(n_bins):
			cg_angle_potentials[angle,i] -= ib_lambda*(ener_temp[angle,i] - atom_angle_potentials[angle,i])
		min_val = np.min(cg_angle_potentials[angle,:])
                cg_angle_potentials[angle,:] -= min_val
		for i in range(n_bins):
			# finite difference force:
			if i==0:
				force[i] = -(cg_angle_potentials[angle,i+1]-cg_angle_potentials[angle,i])/delta_angle
			elif i == n_bins-1:
				force[i] = -(cg_angle_potentials[angle,i]-cg_angle_potentials[angle,i-1])/delta_angle
			else:
				force[i] = -(cg_angle_potentials[angle,i+1]-cg_angle_potentials[angle,i-1])/(2.0*delta_angle)
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,angle_min+(i)*delta_angle, cg_angle_potentials[angle,i],force[i]))
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, ib_out, str(angle+1).strip()))
	out.close()
	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (angle_min+(i+0.5)*delta_angle))
		for angle in range(n_angles):
			temp.write("%12.7f" % (cg_angle_potentials[angle,i]))
		temp.write("\n")
	temp.close()

# MM_ANGLES
#
def UpdateBonds(ib_iter, atom_bond_potentials, cg_bond_hists, cg_bond_potentials, param_out):
	global file_update, kT, ib_lambda, delta_bond

	n_bonds = cg_bond_hists.shape[0]
	n_bins = cg_bond_hists.shape[1]
	ener_temp = np.zeros((n_bonds,n_bins),dtype=float)
	force = np.empty(n_bins,dtype=float)
	coeff_mat = np.empty((n_bins,3),dtype=float)
	x_mat = np.empty(n_bins,dtype=float)
	ib_out = "bond" + str(ib_iter+1) + ".ib"
	ener_out = "bond" + str(ib_iter+1) + ".ener"
	out = open(ib_out,'w')
	for bond in range(n_bonds):
		# create probability density
		cg_bond_hists[bond,:] /= (np.sum(cg_bond_hists[bond,:])*delta_bond)

		for i in range(n_bins):
			x = bond_min+(i+0.5)*delta_bond
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			x_mat[i] = x
			coeff_mat[i,2] = x*x
		# write title to output 
		out.write("\nBOND_%s\n" % (str(bond+1)))
		out.write("N %5d\n\n" % (n_bins))
		# convert to energies
		min_val = 0
		count = 0
		for i in range(n_bins):
			if (cg_bond_hists[bond,i]/coeff_mat[i,2]) > thresh :
				ener_temp[bond,i] = -kT*math.log(cg_bond_hists[bond,i]/coeff_mat[i,2])
			else:
				ener_temp[bond,i] = 0

		min_val = np.amin(ener_temp[bond,:])
#		ener_temp[bond,:] -= min_val
#		cg_bond_potentials[bond,:] -= min_val
		# fit unsampled regions
		i=0
		data_count = 0
		data_flag = "false"
		while i < n_bins:

			if ener_temp[bond,i] != 0 and data_flag == "false":
				start = i
				data_flag = "true"
				stop = i+1
			elif ener_temp[bond,i] != 0:
				stop = i+1
			
			# fit previous section if we have enough data
			if ener_temp[bond,i] == 0 and data_flag == "true" and (stop-start) > 4:
				# add data to missing values
				k, rss, rank, s = np.linalg.lstsq(coeff_mat[start:stop],ener_temp[bond,start:stop])
				for j in range(start):
					ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
					force[j] = -2.0*k[2]*coeff_mat[j,1] - k[1]
				for j in range(stop,n_bins):
					ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
					force[j] = -2.0*k[2]*coeff_mat[j,1] - k[1]

				# compute force
				for j in range(start,stop):
					if j > start and j < stop-1:
						force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j-1])/(2*delta_bond)
					elif j == start:
						force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j])/(delta_bond)
					elif j == stop-1:
						force[j] = -(ener_temp[bond,j]-ener_temp[bond,j-1])/(delta_bond)

				data_flag = "false"
				data_count += 1
				break
			# skip 'island' of data
			elif ener_temp[bond,i] == 0 and data_flag == "true":
				data_flag = "false"

			i += 1
		tck = interpolate.splrep(x_mat,ener_temp[bond,:]-min_val,s=n_bins-math.sqrt(2*n_bins))
		ener_temp[bond,:] = interpolate.splev(x_mat, tck, der=0)
		force = -interpolate.splev(x_mat, tck, der=1)
		min_val = np.min(ener_temp[bond,:])
                ener_temp[bond,:] -= min_val
		for i in range(n_bins):
			cg_bond_potentials[bond,i] -= ib_lambda * (ener_temp[bond,i] - atom_bond_potentials[bond,i])
		min_val = np.min(cg_bond_potentials[bond,:])
                cg_bond_potentials[bond,:] -= min_val
		for i in range(n_bins):
			# finite difference force:
			if i==0:
				force[i] = -(cg_bond_potentials[bond,i+1]-cg_bond_potentials[bond,i])/delta_bond
			elif i == n_bins-1:
				force[i] = -(cg_bond_potentials[bond,i]-cg_bond_potentials[bond,i-1])/delta_bond
			else:
				force[i] = -(cg_bond_potentials[bond,i+1]-cg_bond_potentials[bond,i-1])/(2.0*delta_bond)
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*delta_bond, cg_bond_potentials[bond,i],force[i]))
#			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*delta_bond, f2(coeff_mat[i,1]),force[i]))
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, ib_out, str(bond+1).strip()))
	out.close()

	temp = open(ener_out, 'w')
	for i in range(n_bins):
		temp.write("%12.5f" % (bond_min+(i+0.5)*delta_bond))
		for bond in range(n_bonds):
			temp.write("%12.7f" % (cg_bond_potentials[bond,i]))
		temp.write("\n")
	temp.close()
# MM MM MM

def UpdatePotentials(ib_iter, atom_bond_potentials, cg_bond_hists, cg_bond_potentials, atom_angle_potentials, cg_angle_hists, cg_angle_potentials, atom_dihedral_potentials, cg_dihedral_hists, cg_dihedral_potentials):
	global file_update, kT, ib_lambda, param_out_file, delta_bond, angle_min, delta_angle, dihedral_min, delta_dihedral

	param_out = open(param_out_file, 'w')
	UpdateBonds(ib_iter, atom_bond_potentials, cg_bond_hists, cg_bond_potentials, param_out)
	UpdateAngles(ib_iter, atom_angle_potentials, cg_angle_hists, cg_angle_potentials, param_out)
	UpdateDihedrals(ib_iter, atom_dihedral_potentials, cg_dihedral_hists, cg_dihedral_potentials, param_out)
	param_out.close()

# run CG simulation using lammps
def RunCGSim(ib_iter, sim_in_file, type):
	
	new_sim_file = "cg_sim.iter" + str(ib_iter) + "." + str(type) + ".in"
	new_log_file = "cg_sim.iter" + str(ib_iter) + "." + str(type) + ".log"

	command1 = "sed -e 's/NUM/" + str(ib_iter) + "/g' -e 's/TYPE/" + str(type) + "/g' < " + sim_in_file + " > " +  new_sim_file
	command2 = "time mpirun -np 4 /Users/Ryan/Desktop/lammps-17Nov16/src/lmp_mpi -i " + new_sim_file + " > " + new_log_file	
	
	print command1
	os.system(command1)
	print command2
	os.system(command2)
	print "Done with CG simulation for iteration", ib_iter

#############################################################################################################################################################
####################################################              MAIN PROGRAM               ################################################################
#############################################################################################################################################################

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file

if psf_flag == "true":
	# parse psf file
	ParsePsfFile(top_file)
	print "Number of unique atom types:", n_uniq_atom_types
else:
	# Parse data files containing bond, angle and dihedral data
	ParseDataFiles()


print "Number of bonds:", n_bonds
print "Number of unique bonds:", n_uniq_bonds
print "Number of angles:", n_angles
print "Number of unique angles:", n_uniq_angles
print "Number of dihedrals:", n_dihedrals
print "Number of unique dihedrals:", n_uniq_dihedrals

# declare bond, angle and dihedral histograms
bond_min = 0.0
bond_max = 5.5
delta_bond = 0.005
n_bond_bins  = int((bond_max-bond_min)/delta_bond)
atom_bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
atom_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
cg_bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
cg_bond_potentials = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
angle_min = 0.0
angle_max = 180 # inclusive
delta_angle = 1.0 
n_angle_bins  = int((angle_max-angle_min)/delta_angle) + 1
atom_angle_hists = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
atom_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
cg_angle_hists = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
cg_angle_potentials = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
dihedral_min = -179.0 
dihedral_max = 180.0 # inclusive
delta_dihedral = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/delta_dihedral) + 1
atom_dihedral_hists = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
atom_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
cg_dihedral_hists = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
cg_dihedral_potentials = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)

n_bonds = cg_bond_hists.shape[0]
n_angles = cg_angle_hists.shape[0]
n_dihedrals = cg_dihedral_hists.shape[0]
# read atomistic probability distributions and generate initial (inverse Boltzmann) potentials
ReadAtomHists(atom_bond_hists, atom_bond_potentials, atom_angle_hists, atom_angle_potentials, atom_dihedral_hists, atom_dihedral_potentials)

cg_bond_potentials = atom_bond_potentials
cg_angle_potentials = atom_angle_potentials
cg_dihedral_potentials = atom_dihedral_potentials

# loop through IB iterations for angles
for ib_iter in range(n_iter):
	# run CG simulation
	RunCGSim(ib_iter+1, sim_in_file, "angle")
	
	traj_file = "iter" + str(ib_iter+1) + ".angle.out.dcd"
	
	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_hists,cg_angle_hists,cg_dihedral_hists)

	# update potentials
	if ib_iter == n_iter-1:
		param_out = open("GGGGG.30.cg.dihedral1.lammpspar", 'w')
	else:
		param_out = open("GGGGG.30.cg.angle%s.lammpspar" %(ib_iter+2), 'w')

	for bond in range(n_bonds):
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, "bond0.ib", str(bond+1).strip()))
        UpdateAngles(ib_iter, atom_angle_potentials, cg_angle_hists, cg_angle_potentials, param_out)
	for dihedral in range(n_dihedrals):
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, "dihedral0.ib", str(dihedral+1).strip()))
        param_out.close()

# loop through IB iterations for dihedrals
for ib_iter in range(n_iter):
	# run CG simulation
	RunCGSim(ib_iter+1, sim_in_file, "dihedral")
	
	traj_file = "iter" + str(ib_iter+1) + ".dihedral.out.dcd"
	
	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_hists,cg_angle_hists,cg_dihedral_hists)

	# update potentials
	if ib_iter == n_iter-1:
		param_out = open("GGGGG.30.cg.bond1.lammpspar", 'w')
	else:
		param_out = open("GGGGG.30.cg.dihedral%s.lammpspar" %(ib_iter+2), 'w')
	angle_ib_out = "angle" + str(n_iter) + ".ib"
	for bond in range(n_bonds):
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, "bond0.ib", str(bond+1).strip()))
        for angle in range(n_angles):
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, angle_ib_out, str(angle+1).strip()))
	UpdateDihedrals(ib_iter, atom_dihedral_potentials, cg_dihedral_hists, cg_dihedral_potentials, param_out)
        param_out.close()

# loop through IB iterations for bonded
for ib_iter in range(n_iter):
	# run CG simulation
	RunCGSim(ib_iter+1, sim_in_file, "bond")

	traj_file = "iter" + str(ib_iter+1) + ".bond.out.dcd"

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_hists,cg_angle_hists,cg_dihedral_hists)

	if ib_iter == n_iter-1:
		param_out = open("GGGGG.30.cg.all1.lammpspar", 'w')
	else:
		param_out = open("GGGGG.30.cg.bond%s.lammpspar" %(ib_iter+2), 'w')


        UpdateBonds(ib_iter, atom_bond_potentials, cg_bond_hists, cg_bond_potentials, param_out)
        dihedral_ib_out = "dihedral" + str(n_iter) + ".ib"
	angle_ib_out = "angle" + str(n_iter) + ".ib"
	for angle in range(n_angles):
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, "angle_ib_out", str(angle+1).strip()))
	for dihedral in range(n_dihedrals):
		param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, "dihedral_ib_out.ib", str(dihedral+1).strip()))

	param_out.close()



# loop through IB iterations
for k in range(n_iter):
	
	ib_iter = k + n_iter

	traj_file = "iter" + str(k+1) + ".all.out.dcd"
	
	# run CG simulation
	RunCGSim(k+1, sim_in_file, "all")

	# compute CG probability distributions
	AnalyzeCGTraj(top_file, traj_file,cg_bond_hists,cg_angle_hists,cg_dihedral_hists)

	# update potentials
	param_out = open("GGGGG.30.cg.all%s.lammpspar" %(k+2), 'w')
	UpdatePotentials(ib_iter+1, atom_bond_potentials, cg_bond_hists, cg_bond_potentials, atom_angle_potentials, cg_angle_hists, cg_angle_potentials, atom_dihedral_potentials, cg_dihedral_hists, cg_dihedral_potentials)

