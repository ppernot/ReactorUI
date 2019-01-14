#!/bin/bash

# List results files in MC_Output on all titan_x VMs

# Config parameters
source ./Scripts/cl_config.rc

# Get files from MC_Output
for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    # Get IP of VM to ssh
    VMIP=$(nova list | grep titan_$RUN | awk '{ split($12, v, "="); print v[2]}');
    echo "*** Runs done on titan_${RUN} ***"
    CMD="ls /home/debian/MC_Output/*"
    ssh debian@$VMIP $CMD 2> /dev/null
done

