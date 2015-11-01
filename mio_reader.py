#!/usr/bin/env python
# -*- coding: utf-8 -*- 
import numpy as np
import scipy as sci
import sys, os

def main():
	#change these: atom name and max angular momentum
	atom1 = "N"
	atom2 = "C"
	latom1 = 1
	latom2 = 1	
	opposite = False

	if atom1 == atom2:
		case = "homo"
	else:
		case = "hetero" 
		if latom1 > latom2:
			opposite = True

	skf = open("/Users/at712/plato/Data/TightBinding/mio-1-1/" + atom1 + "-" + atom2 + ".skf", "r")

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
	dx = 		varr[0][0]
	nknots = 	varr[0][1]
	epars = 	varr[1]
	mass =		varr[2]

	# radial distance values for knot points
	rad_range = np.array([n*dx for n in range(nknots)])

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

	# Create numpy array storing all the SKT elements
	skt = np.array(varr[istart:iend])

	# Spline values
	sp_nknots = varr[iend + 1][0]
	sp_cutoff = varr[iend + 1][1]

	# Exponential repulsive potential
	sp_rep = varr[iend + 2]

	# Create numpy array storing all the spline coefficients
	sp_coeff = np.array(varr[iend+3:isp_end])

	print sp_coeff

	def splinefunc(r, sp_coeff, sp_cutoff, sp_rep):
		# Return repulsive potentials if r smaller than first spline knot
		if r < sp_coeff[0][0]:
			return np.exp(-sp_rep[0]*r + sp_rep[1]) + sp_rep[2]

		# Return 0 if r falls outside the cutoff range 
		elif r >= sp_cutoff:
			return 0.0

		elif sp_cutoff > r >= sp_coeff[-2][0]:
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
		return 0

	# SKT knots consistency
	if len(skt) != nknots:
		print "WARNING: number of SKT knots not consistent with number of SKT rows parsed!"
	else:
		print "CHECK: number of SKT knots consistent with number of SKT rows parsed."

	# Spline knots consistency
	if len(sp_coeff) != sp_nknots:
		print "WARNING: number of spline knots not consistent with number of spline rows parsed!"
	else:
		print "CHECK: number of spline knots consistent with number of spline rows parsed."

	# maximum angular momenta for these atoms
	l1 = min(latom1,latom2)
	l2 = max(latom1,latom2)

	# order in which the skt file orbitals are entered 
	def_order = [[0, 0], [0, 1], [0, 2], [1, 1], [1, 1], [1, 2], [1, 2], [1, 2], [2, 2], [2, 2], [2, 2]]
	# order of orbitals appearing in this file
	new_order = [val for val in def_order if val[0] <= l1 and val[1] <= l2]

	print "new order = " + str(new_order)

	# Find empty columns in SKT
	empty_table = (skt == 0)
	empty_columns = np.all(empty_table, 0)

	# Create indices according to columns to copy over
	new_indices = np.array([i for i, val in enumerate(empty_columns) if val == False])

	# Slice out all the 0 elements
	skt = skt[:, new_indices]

	if len(skt[0])%2 != 0:
		print "FATAL WARNING: number of SKT columns not divisible by 2! Number of H and S columns inconsistent!"
		exit(0)
	else:
		print "CHECK: number of SKT columns divisible by 2. Number of H and S columns consistent."

	bdtarr = ["format_3"]

	counter = 0

	# symflag stores the previous orbital #1 #2 id, and is used to check whether
	# the next SKT elements belong to the same orbital with different symmetry
	symflag = 0

	for col in range(len(skt[0]), 0, -2):

		# Only write #1 #2 orbital id if a new interaction is entered
		if new_order[counter] != symflag:
			if opposite:
				new_order[counter] = list(reversed(new_order[counter]))
			orbs = str(new_order[counter][0]) + " " + str(new_order[counter][1])
			bdtarr.append(orbs)

			symflag = new_order[counter]
		
		# Add number of knots
		bdtarr.append(nknots)

		# slice SKT elements related to this radial function
		knotpoints = skt[:, col-2:col]
		if new_order[counter][0] > new_order[counter][1]: #this case means that the orbitals were reversed
			for i in range(nknots):
				knotpoints[i][0] = knotpoints[i][0] * -1 
				knotpoints[i][1] = knotpoints[i][1] * -1 
			sktlist = ['{:>10.2f} {:>17.8e} {:>17.8e}'.format(rad_range[i], knotpoints[i][0], knotpoints[i][1]) for i in range(nknots)]
		else:
			sktlist = ['{:>10.2f} {:>17.8e} {:>17.8e}'.format(rad_range[i], knotpoints[i][0], knotpoints[i][1]) for i in range(nknots)]

		[bdtarr.append(line) for line in sktlist]

		if case == "homo" and new_order[counter][1] > new_order[counter][0]: #in homo case need additional block for reversed symmetry
			new_order[counter] = list(reversed(new_order[counter]))
			orbs = str(new_order[counter][0]) + " " + str(new_order[counter][1])
			bdtarr.append(orbs)
			bdtarr.append(nknots)
			for i in range(nknots):
				knotpoints[i][0] = knotpoints[i][0] * -1 
				knotpoints[i][1] = knotpoints[i][1] * -1 
			sktlist = ['{:>10.2f} {:>17.8e} {:>17.8e}'.format(rad_range[i], knotpoints[i][0], knotpoints[i][1]) for i in range(nknots)]
			[bdtarr.append(line) for line in sktlist]

		counter += 1

	# Add pairpotential
	bdtarr.append(int((sp_cutoff+0.0001)/0.02) + 1)

	pptlist = ['{:>10.2f} {:>17.8e}'.format(rad, splinefunc(rad, sp_coeff, sp_cutoff, sp_rep)) for rad in np.arange(0.0, sp_cutoff+0.0001, 0.02)]
	[bdtarr.append(line) for line in pptlist]

	bdt_file = open("/Users/at712/Desktop/" + atom1 + "_" + atom2 + ".bdt", "w")
	for line in bdtarr:
		bdt_file.write(str(line)+'\n')

	bdt_file.close()
	print "Done."


if __name__ == "__main__":
	# Execute the main code if run as a script.
	main()
