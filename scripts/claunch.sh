#!/bin/bash

# Lauch Virtual Machines titan_x (x=1..CL_SIZE) 
# with image VM_IMAGE and flavor VM_FLAVOR
# defined in cl_config.rc

# Config parameters
source ./Scripts/cl_config.rc

# Launch VMs
for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    echo "Creating VM titan_${RUN}"
    VMUUID=$(nova boot --image "${VM_IMAGE}" --flavor "${VM_FLAVOR}" --key-name "${VM_KEY}" \
        "titan_${RUN}" | awk '/id/ {print $4}' | head -n 1);
    until [[ "$(nova show ${VMUUID} | awk '/status/ {print $4}')" == "ACTIVE" ]]; do
        :
    done
    echo "VM ${RUN} ${VMUUID} is active."
done

# Check status
nova list


