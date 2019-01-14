#!/bin/bash

# Lauch reactor runs with specific ids in given WD

for ((i=$2;i<=$3;i++)) ; do 
  ./Scripts/OneRun_Cloud.sh $1 $i 
done


