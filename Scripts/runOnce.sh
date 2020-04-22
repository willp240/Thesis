#!/bin/bash

cfgFileName=NeighbourClusterConfig
appName=NeighbourClusterApp

RAPTORR_DIR="/vols/build/t2k/wparker2/raptorr/"

ScratchDir="/vols/build/t2k/wparker2/raptorr/"
scriptDir=${ScratchDir}
cfgDir=${ScratchDir}/src/dataproc/tpc/ccd/
logDir=${ScratchDir}

ScriptFileName=${scriptDir}/run${processor}.sh
anaCommand="${RAPTORR_DIR}build/bin/${appName} -c ${cfgDir}/${cfgFileName}.cfg -n 10"
echo "#!/bin/bash" > $ScriptFileName
echo "export RAPTORR_DIR=${RAPTORR_DIR}" >> $ScriptFileName
echo "cat ${cfgDir}/${cfgFileName}.cfg" >> $ScriptFileName

echo "cd ${RAPTORR_DIR}" >> $ScriptFileName
echo "source my_setup.sh" >> $ScriptFileName
echo "source setup_raptorr.sh" >> $ScriptFileName
echo "${anaCommand}" >> $ScriptFileName
chmod 744 $ScriptFileName

stdout="-o ${logDir}/run${processor}.log"
stderr="-e ${logDir}/run${processor}.err"

qsubmit="qsub"
qsubmit="${qsubmit} ${stdout} ${stderr} ${ScriptFileName}"
eval $qsubmit
#	eval $anaCommand
	
echo "done!"
