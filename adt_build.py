#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import numpy as np
import scipy as sci
import sys, os

def main(argv):
	try:
		atom = argv[1]
	except:
		print '\tInput error, expected: ./adt_build.py <atomic species>'
		return 2

	Ry = 2.0
        
        ldict = {"H": 0, "C": 1, "O": 1, "N": 1, "S": 2, "P": 2, "Au": 2, "F": 1, "Cl": 1, "Br": 1, "I": 1}

        if (atom == "Au"):
            rcut = "9.69"
        else: 
            rcut = "5.50"

	skf = open("mio-1-1/" + atom + "-" + atom + ".skf", "r")

	# Import skf file into array
	skfarr = []
	# read over the lines
	for line in skf:
		# buffer the line
		buff = line.replace(",", " ").replace("\t", " ").strip().split()

		for i, entry in enumerate(buff):
			try:
				buff[i] = eval(entry)
			except:
				None
		skfarr.append(buff)

	adtarr = []

	# Add atomic species
	adtarr.append(atom)

	# Add core charge. Assumed to be equal to number of electrons.
	adtarr.append(str(skfarr[1][-1] + skfarr[1][-2] + skfarr[1][-3]))

	# Add number of shells
	#nl = skfarr[0][-1]
	nl = ldict[atom]+1
        adtarr.append(str(nl))

	# Add number of orbitals
	if nl == 1:
		norb = 1
	if nl == 2:
		norb = 4
	if nl == 3:
		norb = 9
	
        adtarr.append(str(norb))

	# Add radius of the orbitals. This seems to be 5.5 bohr for all species in the 
	# mio-1-1 set with a polynomial tail.
	adtarr.append(rcut)

	# Add single site terms in the energy
	# PLACEHOLDER
	adtarr.append("0.0")

	# Add list of quantum numbers for each orbital (n, l, m), then
	for l in range(nl):
		for m in range(2 * l + 1):
			strline = ""

			strline += str(nl) + " " + str(l) + " " + str(m) + " " + str(Ry * skfarr[1][2 - l]) + " "
			
			if l == 0:
				strline += str(skfarr[1][-1]) + "  "
			if l == 1:
				strline += str(skfarr[1][-2]/3.0) + "  "
			if l == 2:
				strline += str(skfarr[1][-3]/5.0) + "  "

			strline += rcut

			adtarr.append(strline)

	adtarr.append("")

	# Add overlap matrix onsite elements
	for i in range(norb):
		strline = ""
		for j in range(norb):
			if i == j:
				strline += "1.0   " 
			else:
				strline += "0.0   "
		adtarr.append(strline)

	adtarr.append("")

	# Add hamilton matrix onsite elements
	for i in range(norb):
		strline = ""
		for j in range(norb):
			if i == j:
				if i < 1:
					strline += str(Ry * skfarr[1][2]) + "   " 
				elif i < 4:
					strline += str(Ry * skfarr[1][1]) + "   "
				else:
					strline += str(Ry * skfarr[1][0]) + "   " 
			else:
				strline += "0.0   "
		adtarr.append(strline)

	adtarr.append("")

	# Add single shell Hubbard U (averaged over all shells) and onsite exchange energy (PLACEHOLDER)
	uvalue = 0

	if nl == 1:
		uvalue = skfarr[1][6]
	if nl == 2: 
		uvalue = 0.5 * (skfarr[1][6] + skfarr[1][5])
	if nl == 3: 
		uvalue = 1.0/3.0 * (skfarr[1][6] + skfarr[1][5] + skfarr[1][4])

	adtarr.append(str(Ry * uvalue) + "  -0.1") 

	for line in adtarr:
		print line


	adt_file = open(atom + ".adt", "w")
	for line in adtarr:
		adt_file.write(str(line)+'\n')

	adt_file.close()
	print "Done."

	return 0

if __name__ == "__main__":
	# Execute the main code if run as a script.
	main(sys.argv)
