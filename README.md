# skt_converter
Python script to convert DFTB+ *.skt tables to PLATO *.bdt tables.

This script yet only supports up to l=1 shells.

How to use:
1. Place your *.skt files into a directory called "mio-1-1" here.
2. Run ./bdt_build.py <atom1> <atom2> for all the combinations of 
   atoms you want to convert into bdt files.

   E.g. for C and H, run bdt_build with C C, H H, C H, and H C

3. If you run bdt_build.py for two different atoms with both l>1,
   eg C and N, you will receive a temporary bdt file '*_het.bdt'.

   These bdt files are not yet ready to use, as the entries for sp
   and ps oribitals are wrong. Run ./spps.py <atom1> <atom2> once,
   and it will create the correct bdt files for <atom1>_<atom2>.bdt
   and <atom2>_<atom1>.bdt both.

4. bdt_build is currently hardcoded to attach smooth polynomial tails
   from r ~= 9.36 bohr to 11 bohr. The range can be changed with 
   setting i and rc to different value. If rc is changed, the orbital
   range in adt_build.py has to be changed to rc/2 to reflec that.

5. Finally, run ./adt_build.py <atom1> to create the atomic data files.

6. In your Data/Tightbinding directory, create a new directory and
   copy your new parameters inside, as well as the model.dat file.

Example for getting hydrocarbons parameters:
./bdt_build.py H H
./bdt_build.py H C
./bdt_build.py C H
./bdt_build.py C C
./adt_build.py C
./adt_build.py H

For N C you would need to also run
./spps.py N C
