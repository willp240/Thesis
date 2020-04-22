#!/bin/bash
MAXJOBS=500
processor=SparkKilled

cfgFileName=${processor}Config
appName=${processor}App

initialRun=111111
initialMask=test1
initialMonth=XX

pathfiles="/vols/dune/data/2018/"
RAPTORR_DIR="/vols/build/t2k/wparker2/raptorr"
count=0

#for i in $(<CERN_beam_test_runs_with_spark_flag_ignored.txt)
#do
for i in {1305005..1355078..1}
do
    for month in 11 12
    do
	for mask in beam gain neutron
	do
	    if [ -e "${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root" ]
	    then
		
		string=$(root -l -b -q './scripts/checkTree.C("'${pathfiles}'/'${month}'/raptorr_raw/raptorr_tpcdata_'${mask}'_R'${i}'.root")')
		if [[ $string == *"header exists"* ]]; then

		    number=$(qstat -u ${USER} | wc -l)
		    while [[ ${number} -ge ${MAXJOBS} ]]; do
			number=$(qstat -u ${USER} | wc -l)
			echo "Sleeping to wait for batch jobs to finish..."
			sleep 50
		    done
		    
		    file="${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root"
		    count=1
		    echo "Doing for ${file}"
		    
                    ScratchDir="${pathfiles}/${month}/SparkFind"
		    scriptDir=${ScratchDir}/scripts
		    cfgDir=${ScratchDir}/cfg
		    logDir=${ScratchDir}/log
		    mkdir -p $scriptDir
		    mkdir -p $cfgDir
		    mkdir -p $logDir
		    
		    sed "s/${initialRun}/${i}/g" ${RAPTORR_DIR}/src/dataproc/tpc/ccd/${cfgFileName}.cfg > ${cfgDir}/${cfgFileName}_${i}.cfg
		    sed -i "s/${initialMask}/${mask}/g" ${cfgDir}/${cfgFileName}_${i}.cfg
		    sed -i "s/${initialMonth}/${month}/g" ${cfgDir}/${cfgFileName}_${i}.cfg
		    
		    ScriptFileName=${scriptDir}/runSparkFind_${i}.sh
		    anaCommand="${RAPTORR_DIR}/build/bin/${appName} -c ${cfgDir}/${cfgFileName}_${i}.cfg "
		    echo "#!/bin/bash" > $ScriptFileName
		    echo "export RAPTORR_DIR=${RAPTORR_DIR}" >> $ScriptFileName
		    echo "cat ${cfgDir}/${cfgFileName}_${i}.cfg" >> $ScriptFileName
		    
		    echo "cd ${RAPTORR_DIR}" >> $ScriptFileName
		    echo "source my_setup.sh" >> $ScriptFileName
		    echo "source setup_raptorr.sh" >> $ScriptFileName
		    echo "${anaCommand}" >> $ScriptFileName
		    chmod 744 $ScriptFileName
		    
		    stdout="-o ${logDir}/runSparkFind_${i}.log"
		    stderr="-e ${logDir}/runSparkFind_${i}.err"
		    
		    qsubmit="qsub"
		    qsubmit="${qsubmit} ${stdout} ${stderr} ${ScriptFileName}"
		    command="${ScriptFileName}"
		    #eval ${ScriptFileName}
		    eval $qsubmit
		else
		    echo "Header doesn't exist in ${pathfiles}/${month}/raptorr_raw/raptorr_tpcdata_${mask}_R${i}.root"
		fi
	    fi	
#	    if (( count == 1 ))
#	    then
#		break 3
#	    fi
	done
    done
done

