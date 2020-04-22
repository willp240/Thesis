#!/bin/bash

echo "starting profile RAM"

if [[ $# -ne 1 && $# -ne 2 ]]
then
  echo "Please supply a PID number corresponding to the job you'd like to profile"
  echo "USAGE: source profileJob.sh <PID> <(optional)output.csv>"
  exit 1
fi

JOBID=$1
OUTPUTNAME=profileJob_${JOBID}.csv
if [[ $# -eq 2 ]]
then
  OUTPUTNAME=$2
fi

STARTTIME=`date +%s`

echo "profiling for PID = ${JOBID}" >> ${OUTPUTNAME}

echo "RAM [KB],VRAM [MB],time since beginning [s]" >> ${OUTPUTNAME}
while ps -p ${JOBID} > /dev/null
do
    RAM=`grep VmRSS /proc/${JOBID}/status`
    RAM=${RAM:6}
    RAM=`echo $RAM | rev | cut -c 4- | rev`
    VRAM=`nvidia-smi | grep ${JOBID}`
    SUBSTRING="ND280_MCMC_2019"
    VRAM=${VRAM#*$SUBSTRING}
    VRAM=`echo $VRAM | rev | cut -c 6- | rev`
    TIME=`date +%s`
    RUNTIME=$((TIME-STARTTIME))
    echo "${RAM},${VRAM},${RUNTIME}" >> ${OUTPUTNAME}
    sleep 5
done

