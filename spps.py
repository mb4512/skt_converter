#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy as np
import scipy as sci
import sys, os

def main():
	
	# Change these: hetero atom names
	atom1 = "N"
	atom2 = "H"

	if atom1 == atom2:
		return 0

	# Open bdt files
	bdt1 = open(atom1 + "_" + atom2 + "_het.bdt", "r")

	bdtarr1 = []

	# Initalise as false to catch missing sp or ps entries
	arr1_sp = False
	arr1_ps = False

	for i, line in enumerate(bdt1):
		bdtarr1.append(line)
		if line == '0 1\n':
			arr1_sp = i
		if line == '1 0\n':
			arr1_ps = i

	bdt1.close()

	bdtarr2 = []
	arr2_sp = False
	arr2_ps = False

	bdt2 = open(atom2 + "_" + atom1 + "_het.bdt", "r")
	for i, line in enumerate(bdt2):
		bdtarr2.append(line)
		if line == '0 1\n':
			arr2_sp = i
		if line == '1 0\n':
			arr2_ps = i

	bdt2.close()

	# Quit if no sp or ps orbitals are found
	if False in [arr1_sp, arr1_ps, arr2_sp, arr2_ps]:
		return 1

	arr1_sp += 2
	arr1_ps += 2
	arr2_sp += 2
	arr2_ps += 2

	# Fetch number of knots
	arr1_sp_nknots = int(bdtarr1[arr1_sp - 1])
	arr1_ps_nknots = int(bdtarr1[arr1_ps - 1])
	arr2_sp_nknots = int(bdtarr2[arr2_sp - 1]) 
	arr2_ps_nknots = int(bdtarr2[arr2_ps - 1])

	# Slice out ps_sigma integrals
	PS1 = bdtarr1[arr1_ps:(arr1_ps + arr1_ps_nknots)]
	PS2 = bdtarr2[arr2_ps:(arr2_ps + arr2_ps_nknots)]

	PS1arr = []
	# Convert to np array
	for line in PS1:
		buff = line.replace(",", " ").replace("\t", " ").strip().split()
		for i, entry in enumerate(buff):
			try:
				buff[i] = eval(entry)
			except:
				None
		PS1arr.append(buff)
	PS1arr = np.array(PS1arr)

	PS2arr = []
	# Convert to np array
	for line in PS2:
		buff = line.replace(",", " ").replace("\t", " ").strip().split()
		for i, entry in enumerate(buff):
			try:
				buff[i] = eval(entry)
			except:
				None
		PS2arr.append(buff)
	PS2arr = np.array(PS2arr)

	# Rebuild string lines and swap signs
	sktlist1 = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in PS1arr]
	sktlist2 = ['{:>10.2f} {:>17.8e} {:>17.8e}\n'.format(row[0], row[1], row[2]) for row in PS2arr]

	# Swap out PS entries
	bdtarr1[arr1_ps:(arr1_ps + arr1_ps_nknots)] = sktlist2
	bdtarr2[arr2_ps:(arr2_ps + arr2_ps_nknots)] = sktlist1

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
	main()
