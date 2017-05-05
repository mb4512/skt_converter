#!/bin/bash

declare -a atoms=('H' 'C' 'O' 'S' 'N' 'P')

rm *.adt > /dev/null
rm *.bdt > /dev/null
rm -r output > /dev/null
mkdir output > /dev/null

echo "Building..."
for i in "${atoms[@]}"
do 
  echo "  $i.adt"
  ./adt_build.py $i | tee "${i}.out" > /dev/null
  for j in "${atoms[@]}"
  do
    echo "    ${i}_${j}.bdt"
    ./bdt_build.py $i $j | tee "${i}_${j}.out" > /dev/null
    echo "    ${j}_${i}.bdt"
    ./bdt_build.py $j $i | tee "${j}_${i}.out" > /dev/null
    ./spps.py $i $j | tee "${i}_${j}_spps.out" > /dev/null
  done
done

mv *.out output/ > /dev/null
mv *het.bdt output/ > /dev/null

echo "Done."
