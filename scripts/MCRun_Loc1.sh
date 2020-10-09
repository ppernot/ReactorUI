#!/bin/bash

if [ ! -z "$3" ]
  then
    cd $3
fi

# Clean results repository
# rm MC_Output/*

for ((i=$1;i<=$1+$2;i++)) 
do 
  Scripts/OneRun_Loc.sh $i $3 
done

echo "==================== All done ! ==================="

