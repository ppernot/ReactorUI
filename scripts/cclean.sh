#!/bin/bash

# Remove all data and run files on titan_x VMs
# (useful if VMs are to be reused)

# Config parameters
source ./Scripts/cl_config.rc

# Get files from MC_Output
for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    # Get IP of VM to ssh
    VMIP=$(nova list | grep titan_$RUN | awk '{ split($12, v, "="); print v[2]}');
    echo "*** Cleaning titan_${RUN} ***"
    CMD="rm -r /home/debian/*"
    ssh debian@$VMIP $CMD 2> /dev/null
done

