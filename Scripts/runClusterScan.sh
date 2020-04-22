#!/bin/bash

#Set name of app we will run, and config we will copy from
cfgFileName=NeighbourClusterConfig
appName=NeighbourClusterApp

#In config set above, you should set these values to be the seed, skirt, and gap. These later get searched for and overwritten with different values for the scan
initialSeed=222222
initialSkirt=444444
initialGap=888888
initialMinSize=111111

#set your raptorr dir, should probably use env var we set in setup_raptorr really
RAPTORR_DIR="/vols/build/t2k/wparker2/raptorr/"

#can use this for running once for tests
count=0

#loop over different values of seed, skirt, gap, and min cluster size
for seed in 40.00 80.00 120.00 160.00
do
    for skirt in 20.00 60.00 80.00 100.00
    do
#	if (( $(echo "${seed} -ge ${skirt}" | bc -l) ));
#	    then
	st=$((`echo "$skirt <= $seed"| bc`))
	if [ $st -eq 1 ]
	then
	    
	    for gap in 0 1 2 5 10 
	    do
		for minsize in 5 10 20 50 
		do
		    #Where you want to save files
		    ScratchDir="/vols/t2k/users/wparker2/RAPTORRFiles/ClusteringOptBlurred/"
		    scriptDir=${ScratchDir}/scripts
		    cfgDir=${ScratchDir}/cfg
		    logDir=${ScratchDir}/log
		    mkdir -p $scriptDir
		    mkdir -p $cfgDir
		    mkdir -p $logDir
		    count=1
		    
		    echo "Doing for Seed ${seed}, Skirt ${skirt}, Gap ${gap}, Min Cluster Size ${minsize}"
		    #Overwrite initial seed skirt and gap in copy of config
		    sed "s/${initialSeed}/${seed}/g" ${RAPTORR_DIR}src/dataproc/tpc/ccd/${cfgFileName}.cfg > ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg 
		    sed -i "s/${initialSkirt}/${skirt}/g" ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg
		    sed -i "s/${initialGap}/${gap}/g" ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg
		    sed -i "s/${initialMinSize}/${minsize}/g" ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg
		    
		    #Make script to submit, sourcing setup scripts and running command
		    ScriptFileName=${scriptDir}/run${processor}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.sh
		    anaCommand="${RAPTORR_DIR}build/bin/${appName} -c ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg -n 1"
       
		    echo "#!/bin/bash" > $ScriptFileName
		    echo "export RAPTORR_DIR=${RAPTORR_DIR}" >> $ScriptFileName
		    echo "cat ${cfgDir}/${cfgFileName}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.cfg" >> $ScriptFileName
		    
		    echo "cd ${RAPTORR_DIR}" >> $ScriptFileName
		    echo "source my_setup.sh" >> $ScriptFileName
		    echo "source setup_raptorr.sh" >> $ScriptFileName
		    echo "${anaCommand}" >> $ScriptFileName
		    chmod 744 $ScriptFileName
		    
		    #Save logs
		    stdout="-o ${logDir}/run${processor}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.log"
		    stderr="-e ${logDir}/run${processor}_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.err"
		    
		    #Submit job!
		    qsubmit="qsub"
		    qsubmit="${qsubmit} ${stdout} ${stderr} ${ScriptFileName}"
		    eval ${ScriptFileName}
		    #eval $qsubmit
		    
		    cp testout.root ${ScratchDir}/Output_Seed${seed}_Skirt${skirt}_Gap${gap}_MinSize${minsize}.root
	    
		    #For debug, prolly ignore
#		    if (( count == 1 ))
#		    then
#		    	break 4
#		    fi
		    
		done
	    done
	fi
    done
done
echo "done!"
