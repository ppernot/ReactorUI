#!/bin/bash

if [ ! -z "$2" ]
  then
    cd $2
fi

cd ./Run

# Get specific data for this run
printf -v NUM4 "%04d" $1

echo ' '
echo '================================================'
echo '================== Run' $NUM4 '===================='
echo '================================================'
echo ' '

cp ../MC_Input/Reactions/run_${NUM4}.csv ./reac_params.dat
rm -f  ./Photo/*.dat
cp ../MC_Input/Photoprocs/${NUM4}_*.dat ./Photo/
cd ./Photo
for f in *.dat; do mv "$f" "$(echo $f | sed s/${NUM4}_// )"; done 
cd -

# Run code
./reactor

# Save results
cp ./fracmol_out.dat ../MC_Output/fracmol_${NUM4}.dat
cp ./rrates.out      ../MC_Output/reacs_rates_${NUM4}.dat
cp ./phrates.out     ../MC_Output/photo_rates_${NUM4}.dat
