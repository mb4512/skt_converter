#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy as np
import scipy as sci
import sys, os

def main(argv):
	try:
		atom1, atom2 = argv[1], argv[2]
	except:
		print '\tInput error, expected: spps.py <atomic species 1> <atomic species 2>'
		return 2

	# Dictionary of shell numbers
	ldict = {"H": 0, "C": 1, "O": 1, "N": 1, "S": 2, "P": 2}

	l1 = ldict[atom1]
	l2 = ldict[atom2]

	# define atom1 to be the one with higher l
	if l1 < l2:
		l1, l2 = l2, l1
		atom1, atom2 = atom2, atom1

	fpath1 = atom1 + "_" + atom2 + "_het.bdt"
	fpath2 = atom2 + "_" + atom1 + "_het.bdt"

	def findentry (btddata):
		divisors = {}
		for i, line in enumerate(btddata):
			if len(line) == 4 and line[1] == " ":
				divisors[line[:-1]] = i
		return divisors

	def toarr (arr):
		# Convert to np array
		nparr = []
		for line in arr:
			buff = line.replace(",", " ").replace("\t", " ").strip().split()
			for i, entry in enumerate(buff):
				try:
					buff[i] = eval(entry)
				except:
					None
			nparr.append(buff)
		nparr = np.array(nparr)
		return nparr


	# Open bdt files
	try:
		bdt1 = open(fpath1, "r")
	except:
		print '\tError: cannot find', fpath1, 'file. Have you called bdt_build for', atom1, atom2,'yet?' 
		return 2

	bdtarr1 = []
	for line in bdt1:
		bdtarr1.append(line)
	bdt1.close()

	try:
		bdt2 = open(fpath2, "r")
	except:
		print '\tError: cannot find', fpath2, 'file. Have you called bdt_build for', atom2, atom1,'yet?' 
		return 2
	
	bdtarr2 = []
	for line in bdt2:
		bdtarr2.append(line)
	bdt2.close()

	orbs1 = findentry(bdtarr1)
	orbs2 = findentry(bdtarr2)

	print orbs1
	print orbs2

	# P-S case
	if l1 == 1 and l2 == 0:
		# Slice out sp_sigma integral
		n_knots = int(bdtarr2[orbs2['0 1'] + 1])
		i_start = int(orbs2['0 1'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 0 orbital
		n_knots2 = int(bdtarr1[orbs1['0 0'] + 1])
		i_start2 = int(orbs1['0 0'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

	# P-P case
	if l1 == 1 and l2 == 1:
		# Slice out sp_sigma integral
		n_knots = int(bdtarr2[orbs2['0 1'] + 1])
		i_start = int(orbs2['0 1'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 1 orbital
		n_knots2 = int(bdtarr1[orbs1['0 1'] + 1])
		i_start2 = int(orbs1['0 1'])

		print n_knots2, i_start2

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		if atom1 == atom2:
			bdtarr2 = bdtarr1

	# D-S case
	if l1 == 2 and l2 == 0:
		# Slice out sp_sigma integral
		n_knots = int(bdtarr2[orbs2['0 1'] + 1])
		i_start = int(orbs2['0 1'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 0 orbital
		n_knots2 = int(bdtarr1[orbs1['0 0'] + 1])
		i_start2 = int(orbs1['0 0'])

		print n_knots2, i_start2

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		# Slice out sd_sigma integral
		n_knots = int(bdtarr2[orbs2['0 2'] + 1])
		i_start = int(orbs2['0 2'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in npSP]

		# Find index after 1 0 orbital
		n_knots2 = int(bdtarr1[orbs1['1 0'] + 1])
		i_start2 = int(orbs1['1 0'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '2 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

	# P-D case
	if l1 == 2 and l2 == 1:
		# Slice out sp_sigma integral
		n_knots = int(bdtarr2[orbs2['0 1'] + 1])
		i_start = int(orbs2['0 1'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 1 orbital
		n_knots2 = int(bdtarr1[orbs1['0 1'] + 1])
		i_start2 = int(orbs1['0 1'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		# Slice out sp_sigma integral
		n_knots = int(bdtarr1[orbs1['0 1'] + 1])
		i_start = int(orbs1['0 1'])
		SP = bdtarr1[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 2 orbital
		n_knots2 = int(bdtarr2[orbs2['0 2'] + 1])
		i_start2 = int(orbs2['0 2'])

		bdtarr2.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr2.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr2.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out sd_sigma integral
		n_knots = int(bdtarr2[orbs2['0 2'] + 1])
		i_start = int(orbs2['0 2'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in npSP]

		# Find index after 1 1 orbital
		n_knots2 = int(bdtarr1[orbs1['1 1'] + 1]) 
		i_start2 = int(orbs1['1 1'])

		n_knots2 +=  int(bdtarr1[i_start2 + n_knots2 + 2])

		bdtarr1.insert(i_start2 + n_knots2 + 3 + 0, '2 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 3 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 3 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out pd_sigma integral
		n_knots = int(bdtarr2[orbs2['1 2'] + 1])
		i_start = int(orbs2['1 2'])

		SP = bdtarr2[i_start+2:i_start + n_knots + 2]
		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]
		stringlist.append(str(n_knots) + '\n')

		SP = bdtarr2[i_start + n_knots + 3:i_start + n_knots + 3 + int(bdtarr2[i_start + n_knots + 2])]
		# Convert to np array
		npSP = toarr(SP)

		[stringlist.append('{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2])) for row in npSP]

		# Find index after 2 0 orbital
		n_knots2 = int(bdtarr1[orbs1['2 0'] + 1]) 
		i_start2 = int(orbs1['2 0'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '2 1\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

	# D-D case
	if l1 == 2 and l2 == 2:
	
		# Slice out sp_sigma integral
		n_knots = int(bdtarr2[orbs2['0 1'] + 1])
		i_start = int(orbs2['0 1'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 2 orbital
		n_knots2 = int(bdtarr1[orbs1['0 2'] + 1])
		i_start2 = int(orbs1['0 2'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		# Slice out sp_sigma integral
		n_knots = int(bdtarr1[orbs1['0 1'] + 1])
		i_start = int(orbs1['0 1'])
		SP = bdtarr1[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]

		# Find index after 0 2 orbital
		n_knots2 = int(bdtarr2[orbs2['0 2'] + 1])
		i_start2 = int(orbs2['0 2'])

		bdtarr2.insert(i_start2 + n_knots2 + 2 + 0, '1 0\n')
		bdtarr2.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr2.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out sd_sigma integral
		n_knots = int(bdtarr2[orbs2['0 2'] + 1])
		i_start = int(orbs2['0 2'])
		SP = bdtarr2[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in npSP]

		# Find index after 1 2 orbital
		n_knots2 = int(bdtarr1[orbs1['1 2'] + 1]) 
		i_start2 = int(orbs1['1 2'])

		n_knots2 +=  int(bdtarr1[i_start2 + n_knots2 + 2])

		bdtarr1.insert(i_start2 + n_knots2 + 3 + 0, '2 0\n')
		bdtarr1.insert(i_start2 + n_knots2 + 3 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 3 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out sd_sigma integral
		n_knots = int(bdtarr1[orbs1['0 2'] + 1])
		i_start = int(orbs1['0 2'])
		SP = bdtarr1[i_start+2:i_start + n_knots + 2]

		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in npSP]

		# Find index after 1 2 orbital
		n_knots2 = int(bdtarr2[orbs2['1 2'] + 1]) 
		i_start2 = int(orbs2['1 2'])

		n_knots2 +=  int(bdtarr2[i_start2 + n_knots2 + 2])

		bdtarr2.insert(i_start2 + n_knots2 + 3 + 0, '2 0\n')
		bdtarr2.insert(i_start2 + n_knots2 + 3 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr2.insert(i_start2 + n_knots2 + 3 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out pd_sigma integral
		n_knots = int(bdtarr2[orbs2['1 2'] + 1])
		i_start = int(orbs2['1 2'])

		SP = bdtarr2[i_start+2:i_start + n_knots + 2]
		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]
		stringlist.append(str(n_knots) + '\n')

		SP = bdtarr2[i_start + n_knots + 3:i_start + n_knots + 3 + int(bdtarr2[i_start + n_knots + 2])]
		# Convert to np array
		npSP = toarr(SP)

		[stringlist.append('{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2])) for row in npSP]

		# Find index after 2 0 orbital
		n_knots2 = int(bdtarr1[orbs1['2 0'] + 1]) 
		i_start2 = int(orbs1['2 0'])

		bdtarr1.insert(i_start2 + n_knots2 + 2 + 0, '2 1\n')
		bdtarr1.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr1.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2

		# Slice out pd_sigma integral
		n_knots = int(bdtarr1[orbs1['1 2'] + 1])
		i_start = int(orbs1['1 2'])

		SP = bdtarr1[i_start+2:i_start + n_knots + 2]
		# Convert to np array
		npSP = toarr(SP)

		# Rebuild string lines and swap signs
		stringlist = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2]) for row in npSP]
		stringlist.append(str(n_knots) + '\n')

		SP = bdtarr1[i_start + n_knots + 3:i_start + n_knots + 3 + int(bdtarr1[i_start + n_knots + 2])]
		# Convert to np array
		npSP = toarr(SP)

		[stringlist.append('{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], -row[1], -row[2])) for row in npSP]

		# Find index after 2 0 orbital
		n_knots2 = int(bdtarr2[orbs2['2 0'] + 1]) 
		i_start2 = int(orbs2['2 0'])

		bdtarr2.insert(i_start2 + n_knots2 + 2 + 0, '2 1\n')
		bdtarr2.insert(i_start2 + n_knots2 + 2 + 1, str(n_knots) + '\n')

		for i,line in enumerate(stringlist):
			bdtarr2.insert(i_start2 + n_knots2 + 2 + 2 + i, line)

		# Rebuild orbital indices
		orbs1 = findentry(bdtarr1)
		orbs2 = findentry(bdtarr2)

		print orbs1
		print orbs2
		
		if atom1 == atom2:
			bdtarr2 = bdtarr1

	bdt_file_1 = open(atom1 + "_" + atom2 + ".bdt", "w")
	bdt_file_2 = open(atom2 + "_" + atom1 + ".bdt", "w")

	for line in bdtarr1:
		bdt_file_1.write(str(line))
	bdt_file_1.close()

	for line in bdtarr2:
		bdt_file_2.write(str(line))
	bdt_file_2.close()

	return 0


if __name__ == "__main__":
	# Execute the main code if run as a script.
	main(sys.argv)
