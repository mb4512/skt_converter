# skt_converter
Python script to convert DFTB+ *.skt tables to PLATO *.bdt tables.

This script yet only supports up to l=2 shells.

How to use:
1. Place your *.skt files into a directory called "mio-1-1" here.

2. Run ./buildall.sh to convert also H,C,O,N,P,S SKT into PLATO format.

3. bdt_build is currently hardcoded to attach smooth polynomial tails
   from r ~= 9.36 bohr to 11 bohr. The range can be changed with 
   setting i and rc to different value. If rc is changed, the orbital
   range in adt_build.py has to be changed to rc/2 to reflec that.

4. In your Data/Tightbinding directory, create a new directory and
   copy your new parameters inside, as well as the model.dat file.
