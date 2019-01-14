#!/bin/bash

# Gather results files in MC_Output on all 
# titan_x VMs into local MC_Output

# Config parameters
source ./Scripts/cl_config.rc

# Get files from MC_Output
for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    # Get IP of VM to ssh
    VMIP=$(nova list | grep titan_$RUN | awk '{ split($12, v, "="); print v[2]}');
    # Copy files 
    echo "*** Result files from titan_${RUN} ***"
    scp debian@$VMIP:/home/debian/MC_Output/* MC_Output 2> /dev/null
done

