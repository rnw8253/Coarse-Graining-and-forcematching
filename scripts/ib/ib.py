import sys
import os
import numpy as np
import math
import MDAnalysis

# read the configuration file and populate the global variables
#def ParseConfigFile(cfg_file):
#        global top_file, traj_file, out_file, file_update
#        f = open(cfg_file)
#        file_update = ""
#        for line in f:
#                # first remove comments                                                                                                                                                                                                                                                                                  
#                if '#' in line:
#                        line, comment = line.split('#',1)
#                if '=' in line:
#                        option, value = line.split('=',1)
#                        option = option.strip()
#                        value = value.strip()
#                        # check value                                                                                                                                                                                                                                                                                    
#                        if option.lower()=='topfile':
#                                top_file = value
#                        elif option.lower()=='trajfile':
#                                traj_file = value
#                        elif option.lower()=='outfile':
#                                out_file = value
#                        elif option.lower()=='fileupdate':
#                                file_update = value
#                        else :
#                                print "Option:", option, " is not recognized"
#        f.close()
#        if file_update != "true" and file_update != "false":
#                file_update = "true"
#                print "file_update being assigned default value of true"
def ParsePsfFile(psf_file):
        global type_list, n_uniq_atom_types, uniq_atom_types, atom_types, uniq_bond_atom_types, n_uniq_bonds, bonds, n_bonds, uniq_angle_atom_types, n_uniq_angles, angles, n_angles,  uniq_dihedral_atom_types, n_uniq_dihedrals, dihedrals, n_dihedrals
        f = open(psf_file)
        atom_flag = bond_flag = angle_flag = dihedral_flag = "false"
        atom = bond = angle = dihedral = 0
        bond_count = 0
        atom_types = []
        bonds = []
        uniq_bond_atom_types = []
        uniq_atom_types = []
        angles = []
        uniq_angle_atom_types = []
        dihedrals = []
	type_list = []	
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
						uniq_atom_types.append(uniq_atom_types[atom2])	
			                        same = "true"
                                                break
                        if same == "false":
                                n_uniq_atom_types += 1
				type_list.append(atom_types[atom])
				uniq_atom_types.append(n_uniq_atom_types)	
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
                                        print n_uniq_bonds+1, uniq_bond_atom_types[n_uniq_bonds]
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
#                                       print n_uniq_angles+1, uniq_angle_atom_types[n_uniq_angles]                                                                                                                                                                                                                                                                                                                                                                                                                                                   
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
#        uniq_atom_types = np.asmatrix(uniq_atom_types,dtype=int)


def computePbcDist2(r1,r2,box,box_2):
        dist2 = 0

        for j in range(0,3):
                temp = r1[j]-r2[j]
                if temp < -box_2[j]:
                        temp += box[j]
                elif temp > box_2[j]:
                        temp -= box[j]
                dist2 += temp*temp

        return dist2;
#
def ComputePbcDih(r1,r2,r3,r4,box):

        r12 = np.empty(3,dtype=float)
        r23 = np.empty(3,dtype=float)
        r34 = np.empty(3,dtype=float)
        # compute two vectors (r1-r2 and r3-r2) making sure they are in same box
        for j in range(0,3):
                temp1 = r1[j]-r2[j]
                temp2 = r2[j]-r3[j]
                temp3 = r3[j]-r4[j]
                if temp1 < -box[j]/2.0:
                        temp1 += box[j]
                elif temp1 > box[j]/2.0:
                        temp1 -= box[j]
                if temp2 < -box[j]/2.0:
                        temp2 += box[j]
                elif temp2 > box[j]/2.0:
                        temp2 -= box[j]
                if temp3 < -box[j]/2.0:
                        temp3 += box[j]
                elif temp3 > box[j]/2.0:
                        temp3 -= box[j]
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

def ComputePbcAng(r1,r2,r3,box):

        r21 = np.empty(3,dtype=float)
        r23 = np.empty(3,dtype=float)

        # compute two vectors (r1-r2 and r3-r2) making sure they are in same box                                               
        for j in range(0,3):
                temp1 = r1[j]-r2[j]
                temp2 = r3[j]-r2[j]
                if temp1 < -box[j]/2.0:
                        temp1 += box[j]
                elif temp1 > box[j]/2.0:
                        temp1 -= box[j]
                if temp2 < -box[j]/2.0:
                        temp2 += box[j]
                elif temp2 > box[j]/2.0:
                        temp2 -= box[j]
                r21[j] = temp1
                r23[j] = temp2

		theta = math.acos(np.dot(r21,r23)/(math.sqrt(np.dot(r21,r21))*math.sqrt(np.dot(r23,r23))))*180.0/3.1415926535
        return theta


	
##############################
#######  Main Script  ########
##############################

# read in command line argument                                                                                                                                                                                                                                                                                          
#cfg_file = sys.argv[1]
top_file = sys.argv[1]
traj_file = sys.argv[2]

#kT = 0.59196
# read cfg file                                                                                                                                                                                                                                                                                                          
#ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file
#print "Output data file:", out_file

# parse psf file                                                                                                                                                                                                                                                                                                         
#ParsePsfFile(top_file)

# Declare universe
u = MDAnalysis.Universe(top_file, traj_file)

ParsePsfFile(top_file)

sel = u.select_atoms("all")
atom_positions = np.empty((sel.n_atoms,3),dtype=float)

n_uniq_bonds = len(u.bonds.types())
n_uniq_angles = len(u.angles.types())
n_uniq_dihedrals = len(u.dihedrals.types())

box = [120,120,120]
box_2 = [60,60,60]
# declare nb dist, bond, angle and dihedral histograms
dist_min = 0.0
dist_max = 55.0
delta_dist = 0.2
n_dist_bins  = int((dist_max-dist_min)/delta_dist)
dist_hists = np.zeros((n_uniq_atom_types,n_uniq_atom_types,n_dist_bins),dtype=float)
dist_counts = np.zeros((n_uniq_atom_types,n_uniq_atom_types),dtype=float)

bond_min = 0.0
bond_max = 5.5
delta_bond = 0.005
n_bond_bins  = int((bond_max-bond_min)/delta_bond)
bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
bond_counts = np.zeros(n_uniq_bonds,dtype=float)

angle_min = 0.0
angle_max = 180 # inclusive
delta_angle = 0.25
n_angle_bins  = int((angle_max-angle_min)/delta_angle) + 1
angle_hists = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
angle_counts = np.zeros(n_uniq_angles,dtype=float)

dihedral_min = -179.0
dihedral_max = 180.0 # inclusive
delta_dihedral = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/delta_dihedral) + 1
dihedral_hists = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
dihedral_counts = np.zeros(n_uniq_dihedrals,dtype=float)

n_atoms = len(sel.atoms)
n_mol = 30
n_res_per_mol = 9
mol_resid_list=[]
for i in range(n_mol):
	temp = [i*9+1,i*9+9]
	mol_resid_list.append(temp)
atom_mol = np.zeros(n_atoms, dtype=int)
for i in range(n_atoms):
	for n in range(n_mol):
		if mol_resid_list[n][0] <= sel.atoms[i].resid <= mol_resid_list[n][1]:
			atom_mol[i] = n
step=0
for ts in u.trajectory:
	print step
 	if step >= 0:
		atom_positions = sel.coordinates()
#		for i in range(n_atoms):
#			for j in range(i,n_atoms):
#				if atom_mol[i] != atom_mol[j]:
#					dist = np.sqrt(computePbcDist2(sel.atoms[i].position,sel.atoms[j].position,box,box_2))
#					dist_bin = int((dist-dist_min)/delta_dist)
#					if dist_min < dist < dist_max:
#						dist_hists[uniq_atom_types[i]-1,uniq_atom_types[j]-1,dist_bin] += 1
#						dist_hists[uniq_atom_types[j]-1,uniq_atom_types[i]-1,dist_bin] += 1
#		
#
#		# Compute Bonds
#		n_bonds = bonds.shape[0]
#		n_bins = bond_hists.shape[1]
#		for bond in range(n_bonds):
#			dist = math.sqrt(computePbcDist2(atom_positions[bonds[bond,0]-1,:],atom_positions[bonds[bond,1]-1,:],box,box_2))
#			dist_bin = int((dist-bond_min)/delta_bond)
#			if dist_bin >= 0 and dist_bin < n_bins:
#				bond_hists[bonds[bond,2],dist_bin] += 1
#		
#
#		# Compute Diheral
#		n_dihedrals = dihedrals.shape[0]
#		n_bins = dihedral_hists.shape[1]
#		for dihedral in range(n_dihedrals):
#			dih = ComputePbcDih(atom_positions[dihedrals[dihedral,0]-1,:],atom_positions[dihedrals[dihedral,1]-1,:],atom_positions[dihedrals[dihedral,2]-1,:],atom_positions[dihedrals[dihedral,3]-1,:],box)
#			dih_bin = int((dih-dihedral_min)/delta_dihedral)
#			if dih_bin >= 0 and dih_bin < n_bins:
#				dihedral_hists[dihedrals[dihedral,4],dih_bin] += 1
#
				
		# Compute angles
		n_angles = angles.shape[0]
		n_bins = angle_hists.shape[1]
		for angle in range(n_angles):
			ang = ComputePbcAng(atom_positions[angles[angle,0]-1,:],atom_positions[angles[angle,1]-1,:],atom_positions[angles[angle,2]-1,:],box)
			ang_bin = int((ang-angle_min)/delta_angle)
			if ang_bin >= 0 and ang_bin < n_bins:
				angle_hists[angles[angle,3],ang_bin] += 1
	step += 1

		
## Write Non-bonded FE
#for i in range(n_uniq_atom_types):
#        for j in range(i,n_uniq_atom_types):
#		out = open('non_bond/%s.%s.nb.hist' %(type_list[i],type_list[j]),'w')
#		for k in range(n_dist_bins):
#			out.write("%10.5f %20.2f\n" %(dist_min+(k+0.5)*delta_dist, dist_hists[i,j,k]))
#		out.close
#
#for bond in range(n_uniq_bonds):	
#	filename = "../bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".hist"
#	out = open(filename, 'w')
#	for i in range(n_bond_bins):
#		out.write("%10.5f %20.2f\n" % (bond_min+(i+0.5)*delta_bond, bond_hists[bond,i]))
#	out.close()

for angle in range(n_uniq_angles):
	filename = "../angs/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".hist"
	out = open(filename, 'w')
	for i in range(n_angle_bins):
		out.write("%10.5f %20.2f\n" % (angle_min+(i+0.5)*delta_angle, angle_hists[angle,i]))
	out.close()
	
	
#for dihedral in range(n_uniq_dihedrals):
#	filename = "../dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".hist"
#	out = open(filename, 'w')
#	for i in range(n_dihedral_bins):
#		out.write("%10.5f %20.2f\n" % (dihedral_min+(i+0.5)*delta_dihedral, dihedral_hists[dihedral,i]))
#	out.close
