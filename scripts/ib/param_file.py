import sys
import os
import numpy as np
import math
import MDAnalysis
from scipy import interpolate

##############################################################
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

##############################################################
def FinalizeNonBonded(dist_hists, dist_min, delta_dist, n_dist_bins, n_nb_types, nb_pairs, nb_pair_num, param_out):
	global kT
	out = open("nonbonded.ib",'w')
	force = np.empty(n_dist_bins,dtype=float)
	coeff_mat = np.empty((n_dist_bins,3),dtype=float)
	coeff_mat2 = np.empty((n_dist_bins,2),dtype=float)
	x_mat = np.empty(n_dist_bins,dtype=float)
	ener_temp = np.zeros((n_nb_types,n_dist_bins),dtype=float)
        log_points = np.zeros((n_nb_types,n_dist_bins),dtype=float)
	# first check to see if we have two copies of the same nonbonded pair
	for nb in range(n_nb_types):
		for i in range(n_dist_bins):
			x = dist_min+(i+0.5)*delta_dist
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			coeff_mat[i,2] = x*x

			coeff_mat2[i,0] = 1.0
			#coeff_mat2[i,0] = 1.0

			coeff_mat2[i,1] = 1/x**(2)

			x_mat[i]=x
		# write title to output
		out.write("\nNB_%s\n" % (str(nb+1)))
		out.write("N %5d\n\n" % (n_dist_bins))

		# Read in histogram files 
                filename = "non_bond/"+nb_pairs[nb][0]+"."+nb_pairs[nb][1]+".nb.hist"
                filename2 = "non_bond/"+nb_pairs[nb][0]+"."+nb_pairs[nb][1]+".nb.hist"
                if os.path.isfile(filename):
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                dist_hists[nb,count] = float(val)
                                count += 1
                        inp.close()
                elif os.path.isfile(filename2):
			filename = filename2
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                dist_hists[nb,count] = float(val)
                                count += 1
                        inp.close()
		else :
			print "Non-bonded pair does not exist:",nb_pairs[nb]
			exit()

		# convert to probability density
		dist_hists[nb,:] /= (np.sum(dist_hists[nb,:])*delta_dist)

		# convert to Free Energies
		min_val = 0
                count = 0
                for i in range(n_dist_bins):
                        if dist_hists[nb,i] > thresh:
				ener_temp[nb,i] = -kT*np.log(dist_hists[nb,i]/coeff_mat[i,2])
                        else:
                                ener_temp[nb,i] = 0
		
		# fit unsampled regions                                                                                                                                                                         
                i=0
                data_count = 0
                data_flag = "false"
                while i < n_dist_bins:
                        if ener_temp[nb,i] != 0 and data_flag == "false":
                                start = i
                                data_flag = "true"
                                stop = i+1
			elif ener_temp[nb,i] != 0:
				stop = i+1
			# fit r**(10) to the first few points for excluded volume
			if data_flag == "true" and data_count == 0:
				init_start = start
#				for j in range(start,start+50):
#					log_points[nb,j] = np.log(ener_temp[nb,j]+10) + 10*np.log((j+0.5)*delta_dist+dist_min)
#				current_points = log_points[nb,start:start+20]
				current_points = ener_temp[nb,start:start+10]
				current_x = coeff_mat2[start:start+10]
				k, rss, rank, s = np.linalg.lstsq(current_x, current_points)
#				k[0] = np.exp(k[0])
				# extrapolate
				for j in range(1,start):
#					ener_temp[nb,j] = k[0]*(1/((j*0.5+dist_min)**(10))) 
					ener_temp[nb,j] = k[0]*coeff_mat2[j,0] + k[1]*coeff_mat2[j,1]
 				# smooth region with cubic spline interpolation                                                                                                                                     
				smooth = 2
#				var1 = int(start-4)
				var1 = 1
				var2 = int(start)
				tck = interpolate.splrep(x_mat[var1:var2],ener_temp[nb,var1:var2],s=smooth)
				ener_temp[nb,var1:var2] = interpolate.splev(x_mat[var1:var2], tck, der=0)
				# set the first bin to be very large 
			       	ener_temp[nb,0] = 1000
                       
			# fit previous section if we have enough data                                                                                                                                           
                        if ener_temp[nb,i] == 0 and data_flag == "true" and (stop-start) > 5:
                                # Fit the first well with a quadratic, smooth and extropalate                                                                                                                                  
                                if data_count == 0:
					init_points = ener_temp[nb,start:start+5]
                                        init_x = coeff_mat[start:start+5]
                                        init_start = start
					previous_points = ener_temp[nb,stop-5:stop]
                                        previous_x = coeff_mat[stop-5:stop]
                                        previous_stop = stop
					previous_start = start


                                else:
                                        # fit to harmonic if big enough gap                                                                                                                                     
                                        if (start-previous_stop) > 6:
                                                current_points = ener_temp[nb,start:start+5]
                                                current_x = coeff_mat[start:start+5]
                                                data_y = np.append(previous_points,current_points)
                                                data_x = np.append(previous_x, current_x, axis=0)
                                                k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
                                                # fill in hole                                                                                                                                                  
                                                for j in range(previous_stop,start):
                                                        ener_temp[nb,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                                # smooth region with cubic spline interpolation                                                                                                                 
                                                smooth = 2
                                                tck = interpolate.splrep(x_mat[previous_stop-3:start+3],ener_temp[nb,previous_stop-3:start+3],s=smooth)
                                                ener_temp[nb,previous_stop-3:start+3] = interpolate.splev(x_mat[previous_stop-3:start+3], tck, der=0)
                                        # otherwise fit to more steep quadratic?                                                                                                                                
                                        else:
                                                current_points = ener_temp[nb,start:start+2]
                                                current_x = coeff_mat[start:start+2]
                                                data_y = np.append(previous_points[-2:],current_points)
                                                data_x = np.append(previous_x[-2:],current_x,axis=0)
                                                k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
                                                # fill in hole                                                                                                                                                  
                                                for j in range(previous_stop,start):
                                                        ener_temp[nb,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                                # smooth region with cubic spline interpolation                                                                                                                 
                                                smooth = 2
                                                tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[nb,previous_stop-2:start+2],s=smooth)
                                                ener_temp[nb,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
                                        # save current data into previous                                                                                                                                       
                                        previous_points = ener_temp[nb,stop-5:stop]
                                        previous_x = coeff_mat[stop-5:stop]
                                        previous_stop = stop
					previous_start = start
                                data_flag = "false"
                                data_count += 1
                        # skip 'island' of data                                                                                                                                                                 
                        elif ener_temp[nb,i] == 0 and data_flag == "true":
                                data_flag = "false"

                        i += 1		
		
		# Shift max between 30-45 angstroms to be the zero
		shift = max(ener_temp[nb,int(30/delta_dist):int(45/delta_dist)])	
		for i in range(n_dist_bins):
			if shift == ener_temp[nb,i]:
				shift_index = i
		for i in range(n_dist_bins):
			if ener_temp[nb,i] != 0:
				ener_temp[nb,i] -= shift
		# Set everything after max to be zero
		ener_temp[nb,shift_index:n_dist_bins] = 0

		# smooth the entire FES
		smooth = 2
		tck = interpolate.splrep(x_mat[1:25],ener_temp[nb,1:25],s=smooth)
		ener_temp[nb,1:25] = interpolate.splev(x_mat[1:25], tck, der=0)
		# smooth the entire FES
		smooth = 0.1
		tck = interpolate.splrep(x_mat[0:n_dist_bins],ener_temp[nb,0:n_dist_bins],s=smooth)
		ener_temp[nb,0:n_dist_bins] = interpolate.splev(x_mat[0:n_dist_bins], tck, der=0)

		# Shift max between 30-45 angstroms to be the zero
		shift = max(ener_temp[nb,int(30/delta_dist):int(45/delta_dist)])	
		for i in range(n_dist_bins):
			if shift == ener_temp[nb,i]:
				shift_index = i
		for i in range(n_dist_bins):
			if ener_temp[nb,i] != 0:
				ener_temp[nb,i] -= shift
		# Set everything after max to be zero
		ener_temp[nb,shift_index:n_dist_bins] = 0

		# compute force using finite difference
		for j in range(n_dist_bins):
			if j > 0 and j < n_dist_bins-1:
				force[j] = -(ener_temp[nb,j+1]-ener_temp[nb,j-1])/(2*delta_dist)
			elif j == 0:
				force[j] = -(ener_temp[nb,j+1]-ener_temp[nb,j])/(delta_dist)
			else:
				force[j] = -(ener_temp[nb,j]-ener_temp[nb,j-1])/(delta_dist)
                

		# Write out nb energies to nbs.ib and param_out
		for i in range(n_dist_bins):
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,dist_min+(i+0.5)*delta_dist, ener_temp[nb,i],force[i]))
		param_out.write("pair_coeff %2d %2d %s NB_%s 45.0\n" %(nb_pair_num[nb][0],nb_pair_num[nb][1], "nonbonded.ib", str(nb+1).strip()))
	out.close()
		
	temp = open("nonbonded.ener", 'w')
	for i in range(n_dist_bins):
		temp.write("%12.5f" % (dist_min+(i+0.5)*delta_dist))
		for nb in range(n_nb_types):
			temp.write(" %15.7f" % (ener_temp[nb,i]-min_val))
		temp.write("\n")
	temp.close()



##############################################################
def FinalizeBond(bond_hists, bond_min, delta_bond, n_bond_bins,n_uniq_bonds, param_out):
	global kT
	out = open("bonds.ib",'w')
	force = np.empty(n_bond_bins,dtype=float)
	coeff_mat = np.empty((n_bond_bins,3),dtype=float)
	x_mat = np.empty(n_bond_bins,dtype=float)
	ener_temp = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
	# Compute bond distances
	for bond in range(n_uniq_bonds):
		for i in range(n_bond_bins):
			x = bond_min+(i+0.5)*delta_bond
			coeff_mat[i,0] = 1.0
			coeff_mat[i,1] = x
			coeff_mat[i,2] = x*x
			x_mat[i] = x
		# write title to output
		out.write("\nBOND_%s\n" % (str(bond+1)))
		out.write("N %5d\n\n" % (n_bond_bins))

		# Read in histogram files
		filename = "../bonds/"+uniq_bond_atom_types[bond][0].strip()+"_"+uniq_bond_atom_types[bond][1].strip()+".hist"
		filename2 = "../bonds/"+uniq_bond_atom_types[bond][1].strip()+"_"+uniq_bond_atom_types[bond][0].strip()+".hist"
		if os.path.isfile(filename):
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				bond_hists[bond,count] = float(val)
				count += 1
			inp.close()
		elif os.path.isfile(filename2):
			filename = filename2
			inp = open(filename,'r')
			count = 0
			for line in inp:
				junk, val = line.split()
				bond_hists[bond,count] = float(val)
				count += 1
			inp.close()
		else:
			print "There is no file for bond:", u.bonds.types()[bond]
			exit()

		# create probability density                                                                                                                                                  
                bond_hists[bond,:] /= (np.sum(bond_hists[bond,:])*delta_bond)

		# Convert to Free Energies                                                                                                                                                        
		min_val = 0
		count = 0
		for i in range(n_bond_bins):
                        if bond_hists[bond,i] > thresh :
                                ener_temp[bond,i] = -kT*math.log(bond_hists[bond,i]/coeff_mat[i,2])
                                if count == 0 or ener_temp[bond,i] < min_val:
                                        min_val = ener_temp[bond,i]
                                count += 1
                        else:
                                ener_temp[bond,i] = 0
		
		# fit unsampled regions                                                                                                                                                                         
                i=0
                data_count = 0
                data_flag = "false"
                while i < n_bond_bins:

                        if ener_temp[bond,i] != 0 and data_flag == "false":
                                start = i
                                data_flag = "true"
                                stop = i+1
				peak_stop = i+1
                        elif ener_temp[bond,i] != 0:
                                stop = i+1

                        # fit previous section if we have enough data                                                                                                                                           
                        if ener_temp[bond,i] == 0 and data_flag == "true" and (stop-start) > 5:
                                # Fit the first well with a quadratic, smooth and extropalate                                                                                                                                  
                                if data_count == 0:
                                        init_points = ener_temp[bond,start:start+5]
                                        init_x = coeff_mat[start:start+5]
                                        init_start = start
					previous_points = ener_temp[bond,stop-5:stop]
                                        previous_x = coeff_mat[stop-5:stop]
                                        previous_stop = stop
					previous_start = start
					
					# fit data
					current_points = ener_temp[bond,start:stop]
					current_x = coeff_mat[start:stop]
					k, rss, rank, s = np.linalg.lstsq(current_x, current_points)
					# extrapolate
					for j in range(start):
						ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
					# smooth region with cubic spline interpolation                                                                                                                                     
					smooth = 2
					tck = interpolate.splrep(x_mat[0:n_bond_bins],ener_temp[bond,0:n_bond_bins],s=smooth)
					ener_temp[bond,0:n_bond_bins] = interpolate.splev(x_mat[0:n_bond_bins], tck, der=0)

                                else:
                                        # fit to harmonic if big enough gap                                                                                                                                     
                                        if (start-previous_stop) > 6:
                                                current_points = ener_temp[bond,start:start+5]
                                                current_x = coeff_mat[start:start+5]
                                                data_y = np.append(previous_points,current_points)
                                                data_x = np.append(previous_x, current_x, axis=0)
                                                k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
                                                # fill in hole                                                                                                                                                  
                                                for j in range(previous_stop,start):
                                                        ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                                # smooth region with cubic spline interpolation                                                                                                                 
                                                smooth = 2
                                                tck = interpolate.splrep(x_mat[previous_stop-3:start+3],ener_temp[bond,previous_stop-3:start+3],s=smooth)
                                                ener_temp[bond,previous_stop-3:start+3] = interpolate.splev(x_mat[previous_stop-3:start+3], tck, der=0)
                                        # otherwise fit to more steep quadratic?                                                                                                                                
                                        else:
                                                current_points = ener_temp[bond,start:start+2]
                                                current_x = coeff_mat[start:start+2]
                                                data_y = np.append(previous_points[-2:],current_points)
                                                data_x = np.append(previous_x[-2:],current_x,axis=0)
                                                k, rss, rank, s = np.linalg.lstsq(data_x, data_y)
                                                # fill in hole                                                                                                                                                  
                                                for j in range(previous_stop,start):
                                                        ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]
                                                # smooth region with cubic spline interpolation                                                                                                                 
                                                smooth = 2
                                                tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[bond,previous_stop-2:start+2],s=smooth)
                                                ener_temp[bond,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
                                        # save current data into previous                                                                                                                                       
                                        previous_points = ener_temp[bond,stop-5:stop]
                                        previous_x = coeff_mat[stop-5:stop]
                                        previous_stop = stop
					previous_start = start
                                data_flag = "false"
                                data_count += 1

			# skip 'island' of data                                                                                                                                                                 
                        elif ener_temp[bond,i] == 0 and data_flag == "true":
                                data_flag = "false"

                        i += 1

		current_points = ener_temp[bond,previous_start:previous_stop]
		current_x = coeff_mat[previous_start:previous_stop]
		k, rss, rank, s = np.linalg.lstsq(current_x, current_points)
		for j in range(previous_stop,n_bond_bins):
			ener_temp[bond,j] = k[0]*coeff_mat[j,0] + k[1]*coeff_mat[j,1] + k[2]*coeff_mat[j,2]

		# smooth region with cubic spline interpolation                                                                                                                                     
		smooth = 8
		tck = interpolate.splrep(x_mat[0:n_bond_bins],ener_temp[bond,0:n_bond_bins],s=smooth)
		ener_temp[bond,0:n_bond_bins] = interpolate.splev(x_mat[0:n_bond_bins], tck, der=0)

		min_val = np.min(ener_temp[bond,:])
      		ener_temp[bond,:] -= min_val


		#force[0:n_bond_bins] = interpolate.splev(x_mat[0:n_bond_bins], tck, der=1)
		# compute force using finite difference                                                                                                                                                                                                              
		for j in range(n_bond_bins):
			if j > 0 and j < n_bond_bins-1:
				force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j-1])/(2*delta_bond)
			elif j == 0:
				force[j] = -(ener_temp[bond,j+1]-ener_temp[bond,j])/(delta_bond)
			else:
				force[j] = -(ener_temp[bond,j]-ener_temp[bond,j-1])/(delta_bond)
                

		# Write out bond energies to bonds.ib and param_out
		for i in range(n_bond_bins):
                        out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,bond_min+(i+0.5)*delta_bond, ener_temp[bond,i],force[i]))
		param_out.write("bond_coeff %2d %s BOND_%s\n" %(bond+1, "bonds.ib", str(bond+1).strip()))
	out.close()

	temp = open("bonds.ener", 'w')
        for i in range(n_bond_bins):
                temp.write("%12.5f" % (bond_min+(i+0.5)*delta_bond))
                for bond in range(n_uniq_bonds):
                        temp.write("%15.7f" % (ener_temp[bond,i]-min_val))
                temp.write("\n")
        temp.close()

#############################################################
def FinalizeAngle(angle_hists, angle_min, delta_angle, n_angle_bins, n_uniq_angles, param_out):
	
	force = np.empty(n_angle_bins,dtype=float)
        coeff_mat = np.empty((n_angle_bins,3),dtype=float)
	x_mat = np.empty((n_angle_bins),dtype=float)
	ener_temp = np.empty((n_uniq_angles,n_angle_bins),dtype=float)

	print n_angle_bins
        out = open("angles.ib",'w')
        # first check to see if we have two copies of same angle
	for angle in range(n_uniq_angles):
                for i in range(n_angle_bins):
                        x = angle_min+(i+0.5)*delta_angle
                        coeff_mat[i,0] = 1.0
                        coeff_mat[i,1] = x
                        coeff_mat[i,2] = x*x
			x_mat[i] = x
		# write title to output 
		out.write("\nANG_%s\n" % (str(angle+1)))
                out.write("N %5d\n\n" % (n_angle_bins))

		# Read in histogram files  
                filename = "../angs/"+uniq_angle_atom_types[angle][0].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][2].strip()+".hist"
                filename2 = "../angs/"+uniq_angle_atom_types[angle][2].strip()+"_"+uniq_angle_atom_types[angle][1].strip()+"_"+uniq_angle_atom_types[angle][0].strip()+".hist"
                if os.path.isfile(filename):
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                angle_hists[angle,count] = float(val)
                                count += 1
                        inp.close()
                elif os.path.isfile(filename2):
                        filename = filename2
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                angle_hists[angle,count] = float(val)
                                count += 1
                        inp.close()
		else:
			print "There is no file for angle:", u.angles.types()[angle]
			exit()

		# convert to probability density
		angle_hists[angle,:] /= (np.sum(angle_hists[angle,:])*delta_angle)
                ener_temp[angle,:] = angle_hists[angle,:]
		# convert to Free Energies
		min_val = 0
                count = 0
                for i in range(n_angle_bins):
                        if angle_hists[angle,i] > thresh:
				#ener_temp[angle,i] = -kT*math.log(angle_hists[angle,i])
				ener_temp[angle,i] = -kT*math.log(angle_hists[angle,i]/math.sin(coeff_mat[i,1]*np.pi/180))
				# angle_potentials[angle,i] = ener_temp[angle,i]				
				# ener_temp[angle,i] = -kT*math.log(ener_temp[angle,i])
                                if count == 0 or ener_temp[angle,i] < min_val:
                                        min_val = ener_temp[angle,i]
                                count += 1
                        else:
                                ener_temp[angle,i] = 0
		
                        # make sure the tails are increasing
                        if i > 0 and ener_temp[angle,i] == 0 and ener_temp[angle,i-1] != 0:
                                for j in range(1,4):
                                        if ener_temp[angle,i-j] < ener_temp[angle,i-1-j]:
                                                ener_temp[angle,i-j] = 0
					#	angle_potentials[angle,i-j] = 0           
		j = 0
                while ener_temp[angle,j] == 0:
                        j += 1
                first_nonzero = j
                for j in range(first_nonzero,first_nonzero+3):
                        if ener_temp[angle,j+1] > ener_temp[angle,j]:
                                ener_temp[angle,j] = 0
			#	angle_potentials[angle,j] = 0
		# fit unsampled regions
		i=0
                data_count = 0
                data_flag = "false"
		while i < n_angle_bins:
                        if ener_temp[angle,i] != 0 and data_flag == "false":
                                start = i
                                data_flag = "true"
                                stop = i+1
                	elif ener_temp[angle,i] != 0:
                                stop = i+1
			#else:
			#	print "false"
		
		
			# fit previous section if we have enough data
			if ener_temp[angle,i] == 0 and data_flag == "true" and (stop-start) > 5:
				# smooth region with cubic spline interpolation
				smooth = 2
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
                                                        force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
						# smooth region with cubic spline interpolation
						smooth = 2
                                                tck = interpolate.splrep(x_mat[previous_stop-2:start+2],ener_temp[angle,previous_stop-2:start+2],s=smooth)
                                                ener_temp[angle,previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=0)
                                                force[previous_stop-2:start+2] = interpolate.splev(x_mat[previous_stop-2:start+2], tck, der=1)                                       
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
                                                        force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
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
		if data_flag == "true" and (n_angle_bins-start) > 5:
			init_points = ener_temp[angle,start:start+5]
			init_x = coeff_mat[start:start+5]
			init_start = start
			previous_stop = n_angle_bins
			previous_x = coeff_mat[n_angle_bins-5:n_angle_bins]
			previous_points = ener_temp[angle,n_angle_bins-5:n_angle_bins]
				
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
		for j in range(previous_stop,n_angle_bins):
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
		if previous_stop < n_angle_bins - 4 :
                        tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[angle,previous_stop-4:previous_stop+4],s=smooth)
                        ener_temp[angle,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
                        force[previous_stop-4:previous_stop+4] = -interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=1)
                else:
                        tck = interpolate.splrep(x_mat[previous_stop-4:n_angle_bins],ener_temp[angle,previous_stop-4:n_angle_bins],s=smooth)
                        ener_temp[angle,previous_stop-4:n_angle_bins] = interpolate.splev(x_mat[previous_stop-4:n_angle_bins], tck, der=0)
                        force[previous_stop-4:n_angle_bins] = -interpolate.splev(x_mat[previous_stop-4:n_angle_bins], tck, der=1)
		
		min_val = np.min(ener_temp[angle,:])
                ener_temp[angle,:] -= min_val
			
		# compute force using finite difference
		for j in range(n_angle_bins):
			if j > 0 and j < n_angle_bins-1:
				force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j-1])/(2*delta_angle)
			elif j == 0:
				force[j] = -(ener_temp[angle,j+1]-ener_temp[angle,j])/(delta_angle)
			else:
				force[j] = -(ener_temp[angle,j]-ener_temp[angle,j-1])/(delta_angle)
                

		for i in range(n_angle_bins):
		#	angle_potentials[angle,i] = ener_temp[angle,i]		
			out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,angle_min+(i)*delta_angle, ener_temp[angle,i],force[i]))
		param_out.write("angle_coeff %2d %s ANG_%s\n" %(angle+1, "angles.ib", str(angle+1).strip()))
	out.close()

	temp = open("angles.ener", 'w')
	for i in range(n_angle_bins):
		temp.write("%12.5f" % (angle_min+(i+0.5)*delta_angle))
		for angle in range(n_uniq_angles):
			temp.write("%12.7f" % (ener_temp[angle,i]))
		temp.write("\n")
	temp.close()

########################################################################
def FinalizeDihedral(dihedral_hists, dihedral_min, delta_dihedral, n_dihedral_bins, n_uniq_dihedrals, param_out):

        force = np.empty(n_dihedral_bins,dtype=float)
        coeff_mat = np.empty((n_dihedral_bins,3),dtype=float)
	ener_temp = np.empty((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
	x_mat = np.empty(n_dihedral_bins,dtype=float)

        out = open("dihedrals.ib",'w')
        # first check to see if we have two copies of same dihedral
	for dihedral in range(n_uniq_dihedrals):
                for i in range(n_dihedral_bins):
                        x = dihedral_min+(i+0.5)*delta_dihedral
                        coeff_mat[i,0] = 1.0
                        coeff_mat[i,1] = x
			x_mat[i] = x
                        coeff_mat[i,2] = x*x
                # write title to output
                out.write("\nDIH_%s\n" % (str(dihedral+1)))
                out.write("N %5d DEGREES\n\n" % (n_dihedral_bins))

		# Read in histogram files 
                filename = "../dihs/"+uniq_dihedral_atom_types[dihedral][0].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][3].strip()+".hist"
                filename2 = "../dihs/"+uniq_dihedral_atom_types[dihedral][3].strip()+"_"+uniq_dihedral_atom_types[dihedral][2].strip()+"_"+uniq_dihedral_atom_types[dihedral][1].strip()+"_"+uniq_dihedral_atom_types[dihedral][0].strip()+".hist"
                if os.path.isfile(filename):
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                dihedral_hists[dihedral,count] = float(val)
                                count += 1
                        inp.close()
                elif os.path.isfile(filename2):
			filename = filename2
                        inp = open(filename,'r')
                        count = 0
                        for line in inp:
                                junk, val = line.split()
                                dihedral_hists[dihedral,count] = float(val)
                                count += 1
                        inp.close()
		else :
			print "Dihedral does not exist:",u.dihedrals.types()[dihedral]
			exit()

		# convert to probability density
		dihedral_hists[dihedral,:] /= (np.sum(dihedral_hists[dihedral,:])*delta_dihedral)
                
		# convert to Free Energies
		min_val = 0
                count = 0
                for i in range(n_dihedral_bins):
                        if dihedral_hists[dihedral,i] > thresh:
				ener_temp[dihedral,i] = -kT*math.log(dihedral_hists[dihedral,i])
                                #dihedral_hists[dihedral,i] = -kT*math.log(dihedral_hists[dihedral,i])
                                if count == 0 or dihedral_hists[dihedral,i] < min_val:
                                        min_val = dihedral_hists[dihedral,i]
                                count += 1
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
		#		dihedral_potentials[dihedral,j] = 0

		# fit unsampled regions
                i=0
                data_count = 0
                data_flag = "false"
                while i < n_dihedral_bins:
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
#                                                        force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
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
#                                                        force[j] = -2*k[2]*coeff_mat[j,1] - k[1]
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
                                        elif j > stop-1:
                                                force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j])/(delta_dihedral)
                                        else:
                                                force[j] = -(ener_temp[dihedral,j]-ener_temp[dihedral,j-1])/(delta_dihedral)

                                data_flag = "false"
                                data_count += 1
                        
			# skip 'island' of data
			elif ener_temp[dihedral,i] == 0 and data_flag == "true":
                                data_flag = "false"

                        i += 1
		
		if data_flag == "true" and (n_dihedral_bins-start) > 5:
                        init_points = ener_temp[dihedral,start:start+5]
                        init_x = coeff_mat[start:start+5]
                        init_start = start
                        previous_stop = n_dihedral_bins
                        previous_x = coeff_mat[n_dihedral_bins-5:n_dihedral_bins]
                        previous_points = ener_temp[dihedral,n_dihedral_bins-5:n_dihedral_bins]
 
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
                for j in range(previous_stop,n_dihedral_bins):
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
		if previous_stop < n_dihedral_bins - 4 :
                        tck = interpolate.splrep(x_mat[previous_stop-4:previous_stop+4],ener_temp[dihedral,previous_stop-4:previous_stop+4],s=smooth)
                        ener_temp[dihedral,previous_stop-4:previous_stop+4] = interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=0)
                        force[previous_stop-4:previous_stop+4] = -interpolate.splev(x_mat[previous_stop-4:previous_stop+4], tck, der=1)
                else:
                        tck = interpolate.splrep(x_mat[previous_stop-4:n_dihedral_bins],ener_temp[dihedral,previous_stop-4:n_dihedral_bins],s=smooth)
                        ener_temp[dihedral,previous_stop-4:n_dihedral_bins] = interpolate.splev(x_mat[previous_stop-4:n_dihedral_bins], tck, der=0)
                        force[previous_stop-4:n_dihedral_bins] = -interpolate.splev(x_mat[previous_stop-4:n_dihedral_bins], tck, der=1)
			
		min_val = np.min(ener_temp[dihedral,:])
                ener_temp[dihedral,:] -= min_val

		# compute force using finite difference
		for j in range(n_dihedral_bins):
			if j > 0 and j < n_dihedral_bins-1:
				force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j-1])/(2*delta_dihedral)
			elif j == 0:
				force[j] = -(ener_temp[dihedral,j+1]-ener_temp[dihedral,j])/(delta_dihedral)
			else:
				force[j] = -(ener_temp[dihedral,j]-ener_temp[dihedral,j-1])/(delta_dihedral)
                
		for i in range(n_dihedral_bins):
		#	dihedral_potentials[dihedral,i] = ener_temp[dihedral, i]
                        out.write("%5d%20.5f%20.5f%20.5f\n" % (i+1,dihedral_min+(i)*delta_dihedral, ener_temp[dihedral,i],force[i]))
                param_out.write("dihedral_coeff %2d %s DIH_%s\n" %(dihedral+1, "dihedrals.ib", str(dihedral+1).strip()))

        out.close()
        temp = open("dihedrals.ener", 'w')
        for i in range(n_dihedral_bins):
                temp.write("%12.5f" % (dihedral_min+(i+0.5)*delta_dihedral))
                for dihedral in range(n_uniq_dihedrals):
                        temp.write("%12.7f" % (ener_temp[dihedral,i]))
                temp.write("\n")
        temp.close()

##############################################################################
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
#                                       print n_uniq_angles+1, uniq_angle_atom_types[n_uniq_angles]                                                                                                                                                                                                                                                                                                                                                                    \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
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
                                                if (atom_types[int(dihedrals[dihedral][0])-1] == atom_types[int(dihedrals[dihedral2][0])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][3])-1]) or (atom_types[int(dihedrals[dihedral]\
[0])-1] == atom_types[int(dihedrals[dihedral2][3])-1] and atom_types[int(dihedrals[dihedral][1])-1] == atom_types[int(dihedrals[dihedral2][2])-1] and atom_types[int(dihedrals[dihedral][2])-1] == atom_types[int(dihedrals[dihedral2][1])-1] and atom_types[int(dihedrals[dihedral][3])-1] == atom_types[int(dihedrals[dihedral2][0])-1]):
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



	
###############################################################################
###############################  Main Script  #################################
###############################################################################

# read in command line argument
#cfg_file = sys.argv[1]
top_file = sys.argv[1]
traj_file = sys.argv[2]

# output parameter file
out_file = 'GGGGG.30.cg.lammpspar'
out2_file = 'GGGGG.30.cg.paircoeff'


# read cfg file
#ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file

# Constants
T = 298 # K
kB = 0.001987 # kcal/mol/K
kT = kB*T # kcal/mol
thresh = 1E-4

# Declare universe
u = MDAnalysis.Universe(top_file, traj_file)
ParsePsfFile(top_file)
print type_list
print uniq_atom_types

# Declare the number of unique atoms, bonds, angles and dihedrals
n_uniq_bonds = len(u.bonds.types())
n_uniq_angles = len(u.angles.types())
n_uniq_dihedrals = len(u.dihedrals.types())
nb_pairs = []
nb_pair_num = []
for i in range(n_uniq_atom_types):
	for j in range(i,n_uniq_atom_types):
		temp = [type_list[i],type_list[j]]
		temp2 = [i+1,j+1]
		nb_pairs.append(temp)
		nb_pair_num.append(temp2)
n_nb_types = len(nb_pairs)

# Box Size
box = [120,120,120]
box_2 = [60,60,60]

# declare non-bonded, dist, bond, angle and dihedral histograms
dist_min = 0.0
dist_max = 55.0
delta_dist = 0.2
n_dist_bins  = int((dist_max-dist_min)/delta_dist)
dist_hists = np.zeros((n_nb_types,n_dist_bins),dtype=float)
#dist_counts = np.zeros((n_nb_types),dtype=float)

bond_min = 0.0
bond_max = 5.5
delta_bond = 0.005
n_bond_bins  = int((bond_max-bond_min)/delta_bond)
bond_hists = np.zeros((n_uniq_bonds,n_bond_bins),dtype=float)
#bond_counts = np.zeros(n_uniq_bonds,dtype=float)

angle_min = 0.0
angle_max = 180 # inclusive
delta_angle = 1.0
n_angle_bins  = int((angle_max-angle_min)/delta_angle) + 1
angle_hists = np.zeros((n_uniq_angles,n_angle_bins),dtype=float)
#angle_counts = np.zeros(n_uniq_angles,dtype=float)

dihedral_min = -179.0
dihedral_max = 180.0 # inclusive
delta_dihedral = 1.0
n_dihedral_bins  = int((dihedral_max-dihedral_min)/delta_dihedral) + 1
dihedral_hists = np.zeros((n_uniq_dihedrals,n_dihedral_bins),dtype=float)
#dihedral_counts = np.zeros(n_uniq_dihedrals,dtype=float)

param_out =open(out_file,'w')
param2_out=open(out2_file,'w')
# Finish Bond Probability Distributions
FinalizeBond(bond_hists, bond_min, delta_bond, n_bond_bins, n_uniq_bonds, param_out)
# Finish Angle Probibility Distributions
FinalizeAngle(angle_hists, angle_min, delta_angle, n_angle_bins, n_uniq_angles, param_out)
# Finish Dihedral Probibility Distributions
FinalizeDihedral(dihedral_hists, dihedral_min, delta_dihedral, n_dihedral_bins, n_uniq_dihedrals, param_out)
# Finish Nonbonded Probability Distributions
#FinalizeNonBonded(dist_hists, dist_min, delta_dist, n_dist_bins, n_nb_types, nb_pairs, nb_pair_num, param2_out)

param_out.close()
