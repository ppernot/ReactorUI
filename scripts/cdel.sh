#!/bin/bash

# Delete all titan_x VMs

# Config parameters
source ./Scripts/cl_config.rc

for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    nova delete titan_${RUN}
done

# Check
nova list

