#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy as np
import scipy as sci
import sys, os

from polytail import Tail

def main(argv):
	try:
		atom1, atom2 = argv[1], argv[2]
	except:
		print '\tInput error, expected: bdt_build.py <atomic species 1> <atomic species 2>'
		return 2

	# Units conversion
	Ry = 2.0;

	# Dictionary of shell numbers
	ldict = {"H": 0, "C": 1, "O": 1, "N": 1, "S": 2, "P": 2}

	latom1 = ldict[atom1]
	latom2 = ldict[atom2]

	opposite = False

	if atom1 == atom2:
		case = "homo"
	else:

		case = "hetero"

		if latom1 > latom2:
			opposite = True

	skf = open("mio-1-1/" + atom1 + "-" + atom2 + ".skf", "r")

	varr = []
	# read over the lines
	for line in skf:
		# buffer the line
		buff = line.replace(",", " ").replace("\t", " ").strip().split()

		for i, entry in enumerate(buff):
			try:
				buff[i] = eval(entry)
			except:
				None

		varr.append(buff)


	# SKT values
	dx 		= 	varr[0][0]
	nknots  = 	varr[0][1]
	epars   = 	varr[1]
	mass    =	varr[2]

	# find start of tables
	istart = 0

	iend = len(varr)

	for i in range(3, len(varr)):
		# find first occurence after line 2 in which more than 1 column appears (sk tables)
		if istart == 0 and len(varr[i]) > 1:
			istart = i

		# find 'Spline' (end of sk tables)
		if varr[i] == ['Spline']:
			iend = i

		# find '<Documentation>' (end of spline tables)
		if varr[i] == ['<Documentation>']:
			isp_end = i

	# Fetch r offset
	#offset = (istart - 4) * dx
	#print "noffset: ", offset

	# Define r offsset
	offset = 0.38
	print "noffset ", 0.38

	# radial distance values for knot points
	rad_range = np.array([(n+1)*dx for n in range(nknots)]) + offset

	# Create numpy array storing all the SKT elements
	# DEBUG DEBUG DEBUG: the -1
	skt = np.array(varr[istart:iend-1])

	# Spline values
	sp_nknots = varr[iend + 1][0]
	sp_cutoff = varr[iend + 1][1]

	# Exponential repulsive potential
	sp_rep = varr[iend + 2]

	# Create numpy array storing all the spline coefficients
	sp_coeff = np.array(varr[iend+3:isp_end])


	def splinefunc(r, sp_coeff, sp_cutoff, sp_rep):
		# Return repulsive potentials if r smaller than first spline knot
		if r < sp_coeff[0][0]:
			return np.exp(-sp_rep[0]*r + sp_rep[1]) + sp_rep[2]

		# Return 0 if r falls outside the cutoff range 
		elif r >= sp_cutoff:
			return 0.0

		elif sp_cutoff > r >= sp_coeff[-1][0]:
			tailco = sp_coeff[-1]
			dr = r - tailco[0]
			# Lazy polynomial :)
			return tailco[2] + tailco[3]*dr + tailco[4]*np.power(dr, 2) + tailco[5]*np.power(dr, 3) + tailco[6]*np.power(dr, 4) + tailco[7]*np.power(dr, 5)

		# Otherwise return the spline
		else:
			# Loop through all coefficients except the tal
			for coeff in sp_coeff[:-1]:
				if coeff[1] > r >= coeff[0]:
					dr = r - coeff[0]
					return coeff[2] + coeff[3]*dr + coeff[4]*np.power(dr, 2) + coeff[5]*np.power(dr, 3)

		# If everything fails, return 0 and throw an exception
		print "WARNING: Spline value {:f} fell through the spline function!" % r
		return 2

	def makeknots (col, new_order, skt, dx, i, rad_range, n):
		h_index = col
		s_index = col - len(skt[0])/2

		print "s, h: ", s_index, h_index
		# Fetch SKT elements related to this radial function. First column: S, second: H
		knotpoints = np.array([skt[:, h_index], skt[:, s_index]]).T

		# Initialise tail function going from ~9.36 to 11.0 bohr
		i = 450
		r0 = rad_range[i]
		rc = 11.0

		tailH = Tail(knotpoints[:, 0], i, dx, r0, rc)
		tailS = Tail(knotpoints[:, 1], i, dx, r0, rc)

		# Extend rad_range to new lengths
		tail_range = r0 + np.array([n*dx for n in range(int((rc-r0)/0.02) + 2)])
		new_rrange = np.append(rad_range[:i], tail_range)

		# Calculate tail elements
		tailpointsH = np.array([tailH.tail(x) for x in tail_range])
		tailpointsS = np.array([tailS.tail(x) for x in tail_range])

		# Append the tail to knotpoints 
		tailarray = np.array([tailpointsH, tailpointsS]).T
		knotpoints = np.concatenate((knotpoints[:i], tailarray), axis = 0)

		return (len(knotpoints), new_rrange, knotpoints)


	# SKT knots consistency
	if len(skt) != nknots:
		print "WARNING: number of SKT knots not consistent with number of SKT rows parsed!"
		print len(skt), nknots
		nknots = len(skt)
	else:
		print "CHECK: number of SKT knots consistent with number of SKT rows parsed."

	# Spline knots consistency
	if len(sp_coeff) != sp_nknots:
		print "WARNING: number of spline knots not consistent with number of spline rows parsed!"
	else:
		print "CHECK: number of spline knots consistent with number of spline rows parsed."


	# maximum angular momenta for these atoms
	l1 = latom1
	l2 = latom2

	# order in which the skt file orbitals are entered 
	def_order = [[0, 0], [0, 1], [0, 2], [1, 1], [1, 1], [1, 2], [1, 2], [2, 2], [2, 2], [2, 2]]
	
	# order of orbitals appearing in this file
	new_order = [val for val in def_order if (val[0] <= l1 and val[1] <= l2)]
	print "new_order", new_order
	print "Number of SKT columns to import: " + str(2 * len(new_order))

	# Find empty columns in SKT
	empty_table = (skt == 0)
	empty_columns = np.all(empty_table, 0)

	# Create indices according to columns to copy over
	if empty_columns.shape != ():
		new_indices = np.array([i for i, val in enumerate(empty_columns) if val == False])
	else:
		new_indices = range(len(def_order))

	# Slice out all the 0 elements
	skt = skt[:, new_indices]

	print "Shells to import: ", new_order
	print "\nFirst row of skt table: "

	# The different symmetric cases of integrals (example: pp_sigma, pp_pi)	are entered in reverse
	# order in the skt file compared to plato. Hence these entried need to be swapped in the skt array.
	# Create reverse new_order array
	rev_new_order = new_order[::-1]

	dd_dup = [i for i, x in enumerate(rev_new_order) if x == [2, 2]] # dd duplicates 
	pd_dup = [i for i, x in enumerate(rev_new_order) if x == [1, 2]] # pd duplicates
	pp_dup = [i for i, x in enumerate(rev_new_order) if x == [1, 1]] # pp duplicates

	if len(pp_dup) > 0:
		pp_dup = pp_dup[0]

		# Get number of pp shells 
		pp_len = def_order.count([1, 1])

		# Copy skt table into buffer
		new_skt = skt[:]

		# Reverse order of pp shells in skt

		# Reverse hamiltonian entries
		skt[:, pp_dup:pp_dup+pp_len] = np.fliplr(skt[:, pp_dup:pp_dup+pp_len])

		# Reverse overlap entries
		s1 = pp_dup + len(new_order)
		s2 = pp_dup + pp_len + len(new_order)

		print skt[0, s1:s2]
		skt[:, s1:s2] = np.fliplr(skt[:, s1:s2])
		print skt[0, s1:s2]

	if len(pd_dup) > 0:
		pd_dup = pd_dup[0]

		# Get number of pd shells 
		pd_len = def_order.count([1, 2])

		# Reverse hamiltonian entries
		skt[:, pd_dup:pd_dup+pd_len] = np.fliplr(skt[:, pd_dup:pd_dup+pd_len])

		# Reverse overlap entries
		s1 = pd_dup + len(new_order)
		s2 = pd_dup + pd_len + len(new_order)
	
		print skt[0, s1:s2]
		skt[:, s1:s2] = np.fliplr(skt[:, s1:s2])
		print skt[0, s1:s2]

	if len(dd_dup) > 0:
		dd_dup = dd_dup[0]

		# Get number of dd shells 
		dd_len = def_order.count([2, 2])

		# Reverse hamiltonian entries
		skt[:, dd_dup:dd_dup+dd_len] = np.fliplr(skt[:, dd_dup:dd_dup+dd_len])

		# Reverse overlap entries
		s1 = dd_dup + len(new_order)
		s2 = dd_dup + dd_len + len(new_order)

		print skt[0, s1:s2]
		skt[:, s1:s2] = np.fliplr(skt[:, s1:s2])
		print skt[0, s1:s2]

	print "\nFirst row of reversed skt table: "
	print skt[0, :]

	# Check if number of skt columns agree with number to be imported by new_order variable
	if len(skt[0])%2 != 0 and len(skt[0]) != 2*len(new_order):
		print "\nFATAL WARNING: number of SKT columns not divisible by 2! Number of H and S columns inconsistent!"
		exit(0)
	else:
		print "\nCHECK: number of SKT columns divisible by 2. Number of H and S columns consistent."

	bdtarr = ["format_3"]
	counter = 0

	# symflag stores the previous orbital #1 #2 id, and is used to check whether
	# the next SKT elements belong to the same orbital with different symmetry
	symflag = 0

	# loop through the integrals in reverse, starting with ss
	for col in range(len(new_order)):

		# Only write #1 #2 orbital identifier if a new interaction is entered
		if new_order[counter] != symflag:
			orbs = str(new_order[counter][0]) + " " + str(new_order[counter][1])
			bdtarr.append(orbs)
			symflag = new_order[counter]

		(nknots, radrange, knotpoints) = makeknots (len(skt[1])-1 - col, new_order, skt, dx, 450, rad_range, n)
		print col, len(skt[1])-1 - col, knotpoints[0], new_order[counter]
		bdtarr.append(nknots)

		sktlist = ['{:>10.2f} {:>17.8e} {:>17.8e}'.format(radrange[i], knotpoints[i][0], Ry * knotpoints[i][1]) for i in range(nknots)]

		[bdtarr.append(line) for line in sktlist]
		counter += 1

	# Add pairpotential
	bdtarr.append(int((sp_cutoff + 0.0001)/0.02) + 1)

	pptlist = ['{:>10.2f} {:>17.8e}'.format(rad, Ry * splinefunc(rad, sp_coeff, sp_cutoff, sp_rep)) for rad in np.arange(0.0, sp_cutoff+0.0001, 0.02)]
	[bdtarr.append(line) for line in pptlist]


	bdt_file = open(atom1 + "_" + atom2 + "_het.bdt", "w")

	for line in bdtarr:
		bdt_file.write(str(line)+'\n')

	bdt_file.close()
	print "Done."

if __name__ == "__main__":
	# Execute the main code if run as a script.
	main(sys.argv)



