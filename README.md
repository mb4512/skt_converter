# skt_converter
Python script to convert DFTB+ *.skt tables to PLATO *.bdt tables.

This script yet only supports up to l=2 shells.

How to use:
1. Place your *.skt files into a directory called "mio-1-1" in the
   directory of this script. The skt directory can be called
   differently, but then you will need to change the paths in 
   bdt_build.py and adt_build.py

2. Run ./buildall.sh to convert all H,C,N,O,S,P SKT into PLATO format.

3. bdt_build is currently hardcoded to attach smooth polynomial tails
   from r ~= 9.36 bohr to 11 bohr (with some exceptions, see Au). 
   The range can be changed by setting i and rc to different values. 
   If rc is changed, the orbital range in adt_build.py has to be 
   changed to rc/2 to reflect that.

4. In your Data/Tightbinding directory, create a new directory and
   copy your new parameters inside, as well as a new model.dat file.
