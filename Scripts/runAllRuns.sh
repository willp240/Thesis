#!/bin/bash
processor=GaussianBlurring
inputprocessed=SparkKilled

cfgFileName=${processor}Config
appName=${processor}App

initialRun=111111
initialMask=test1
initialMonth=XX

pathfiles="/vols/t2k/users/wparker2/RAPTORRFiles"
pathoutfiles="/vols/dune/data/2018/"
RAPTORR_DIR="/vols/build/t2k/wparker2/raptorr"
count=0
MAXJOBS=400

for i in 1323061 # $(<scripts/runlists/CERN_beam_test_runs_with_spark_flag_ignored.txt)
do
#for i in $(<scripts/runlists/allRuns11_12.txt)
#do
    for month in 08 09 10 11 12
    do
	for mask in beam gain neutron
	do

	    number=$(qstat -u ${USER} | wc -l)
	    while [[ ${number} -ge ${MAXJOBS} ]]; do
		number=$(qstat -u ${USER} | wc -l)
		echo "Sleeping to wait for batch jobs to finish..."
		sleep 50
            done

	    if [ -e "${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root" ]
#	    if [ -e "${pathfiles}/${month}/${inputprocessed}NoFirstBias/${inputprocessed}R${i}.root" ]
	    then

		if [ -e "${pathoutfiles}/${month}/${processor}/${processor}R${i}.root" ]
		then
		    echo "Already done for Run $i"
		    
		else

#		    string=$(root -l -b -q './scripts/checkTree.C("'${pathfiles}'/'${month}'/raptorr_raw/raptorr_tpcdata_'${mask}'_R'${i}'.root")')
		    string=$(root -l -b -q './scripts/checkTree.C("'${pathfiles}'/'${month}'/'${inputprocessed}'NoFirstBias/'${inputprocessed}'R'${i}'.root")')
		    if [[ $string == *"header exists"* ]]; then
			
			#file="${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root"
			file="${pathfiles}/${month}/${inputprocessed}/${inputprocessed}R${i}.root"
			count=1
			
			ScratchDir="${pathoutfiles}/${month}/${processor}"
			scriptDir=${ScratchDir}/scripts
			cfgDir=${ScratchDir}/cfg
			logDir=${ScratchDir}/log
			mkdir -p $scriptDir
			mkdir -p $cfgDir
			mkdir -p $logDir
			
			sed "s/${initialRun}/${i}/g" ${RAPTORR_DIR}/src/dataproc/tpc/ccd/${cfgFileName}.cfg > ${cfgDir}/${cfgFileName}_${i}.cfg
			sed -i "s/${initialMask}/${mask}/g" ${cfgDir}/${cfgFileName}_${i}.cfg
			sed -i "s/${initialMonth}/${month}/g" ${cfgDir}/${cfgFileName}_${i}.cfg
			
			ScriptFileName=${scriptDir}/run${processor}_${i}.sh
			anaCommand="${RAPTORR_DIR}/build/bin/${appName} -c ${cfgDir}/${cfgFileName}_${i}.cfg "
			echo "#!/bin/bash" > $ScriptFileName
			echo "export RAPTORR_DIR=${RAPTORR_DIR}" >> $ScriptFileName
			echo "cat ${cfgDir}/${cfgFileName}_${i}.cfg" >> $ScriptFileName
			
			echo "cd ${RAPTORR_DIR}" >> $ScriptFileName
			echo "source my_setup.sh" >> $ScriptFileName
			echo "source setup_raptorr.sh" >> $ScriptFileName
			echo "${anaCommand}" >> $ScriptFileName
			chmod 744 $ScriptFileName
			
			stdout="-o ${logDir}/run${processor}_${i}.log"
			stderr="-e ${logDir}/run${processor}_${i}.err"
			
			qsubmit="qsub"
			qsubmit="${qsubmit} -l h_vmem=20G ${stdout} ${stderr} ${ScriptFileName}"
			command="${ScriptFileName}"
			#eval ${ScriptFileName}
			eval $qsubmit
			
		    else
#			echo "Header doesn't exist in ${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root"
			echo "Header doesn't exist in ${pathfiles}/${month}/${inputprocessed}NoFirstBias/${inputprocessed}R${i}.root"
		    fi
		fi
	    fi
#	    if (( count == 1 ))
#	    then
#		break 3
#	    fi
	done
    done
done

