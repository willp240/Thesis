#!/bin/bash

for run in run2a run2w run3 run4a run4w run5 run6 run7 run8a run8w run9; do
#for run in run2a ; do    
    input_dir=/vols/t2k/users/wparker2/Throws2020_v2r39/${run}/root
    output_dir=/vols/t2k/users/wparker2/NIWGReWeights2020_v2r39/${run}
    
    mkdir -p ${output_dir}/root
    mkdir -p ${output_dir}/log
    
    # Setup up some basics
    date=$(date +%m_%d)
    timereq="0:29:59"
    MAXJOBS=500
    qsubopt=" -cwd -V -l h_rt=${timereq}"
    
    counter=0
    total=$(find ${input_dir} -name "*.root" -type f | wc -l)
    
    # Run over some files
    for file in $(find ${input_dir} -type f); do
	# Increment the counter
	counter=$((counter+1))
	
	echo "Running on file $i which is ${counter}/${total} of ${run}"
	# Only query if we've submitted 2000 jobs already
	if [[ ${counter} -ge ${MAXJOBS} ]]; then
	    exit 1
	    number=$(qstat -u ${USER} | wc -l)
	    while [[ ${number} -ge ${MAXJOBS} ]]; do
		number=$(qstat -u ${USER} | wc -l)
		echo "Sleeping 50s to wait for batch jobs to finish..."
		sleep 500
	    done
	fi
	
	# Get the base name
	short=$(basename ${file})
	
	output="${output_dir}/root/${short%%_12_14_Throws.root}_${date}_NIWGWeights.root"
	
	log="${output_dir}/log/${short%%_12_14_Throws.root}_${date}_NIWGWeights.log"
	logerr="${output_dir}/log/${short%%_12_14_Throws.root}_${date}_NIWGWeights.err"
	
	cmd="qsub ${qsubopt} -o ${log} -e ${logerr} ./RunThemWeights_wrapper.sh ${file} ${output}"
	#echo $cmd
	eval $cmd
	
    done
done
