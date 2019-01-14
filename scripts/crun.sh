#!/bin/bash

# Copy datasets and run MC loop for 
# (sequential) prescribed indices 
# on each titan_x VM

# Config parameters
source ./Scripts/cl_config.rc

# Clean project's results repository
rm MC_Output/* 2>/dev/null

# Distribute MC runs on VMs
RUNS_BY_VM=$((MC_RUNS/CL_SIZE))
RUNS_BY_CORE=$((RUNS_BY_VM/NB_CORE))
MC_START=0
for ((RUN=1;RUN<=$CL_SIZE;RUN++)) ; do
    # Get IP of VM to ssh
    VMIP=$(nova list | grep titan_$RUN | awk '{ split($12, v, "="); print v[2]}');
    #ssh-keygen -f "/home/pernot/.ssh/known_hosts" -R $VMIP

    # Copy data and scripts files
    echo "Copying database and scripts on titan_${RUN}" 
    scp -r MC_*    debian@$VMIP:/home/debian 2>/dev/null >/dev/null
    scp -r Scripts debian@$VMIP:/home/debian 2>/dev/null >/dev/null

    echo "Launching Reactor on titan_${RUN}"
    for ((CORE=1;CORE<=$NB_CORE;CORE++)) ; do

        # Config Working directory
        scp -r Run  debian@$VMIP:/home/debian/Run_$CORE 2>/dev/null >/dev/null
        
        # Run
        MC_STOP=$((MC_START+RUNS_BY_CORE-1))
        CMD="/home/debian/Scripts/MCRun_Cloud.sh $CORE $MC_START $MC_STOP"
        ssh debian@$VMIP $CMD  </dev/null >/dev/null 2>&1 &

        # Next
        MC_START=$((MC_STOP+1))
    done
done

