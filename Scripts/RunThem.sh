#!/bin/bash

input_dir=/vols/t2k/users/cvw09/data/procs/merged/run2a
output_dir=/vols/t2k/users/cvw09/data/procs/procs/cov/run2a

mkdir -p ${output_dir}/root
mkdir -p ${output_dir}/log

# Setup up some basics
date=$(date +%m_%d)
timereq="0:29:59"
MAXJOBS=1800
qsubopt="-q hep.q -cwd -V -m a -M c.wret14@imperial.ac.uk -l h_rt=${timereq}"

counter=0
total=$(find ${input_dir} -name "*.root" -type f | wc -l)

# Run over some files
#for file in $(find ${input_dir} -name "oa_nt_beam_90300109-0054_wpkgkmin6ayv_anal_000_magnet201011airrun3-bsdv01_2.root" -type f); do
for file in $(find ${input_dir} -type f); do
  # Increment the counter
  counter=$((counter+1))

  echo "Running on file $i which is ${counter}/${total}"
  # Only query if we've submitted 2000 jobs already
  if [[ ${counter} -ge ${MAXJOBS} ]]; then
    number=$(qstat -u ${USER} | wc -l)
    while [[ ${number} -ge ${MAXJOBS} ]]; do
      number=$(qstat -u ${USER} | wc -l)
      echo "Sleeping 50s to wait for batch jobs to finish..."
      sleep 50
    done
  fi

  # Get the base name
  short=$(basename ${file})

  output="${output_dir}/root/${short%%.root}_${date}_FlatProc.root"

  log="${output_dir}/log/${short%%.root}_${date}_FlatProc.log"
  logerr="${output_dir}/log/${short%%.root}_${date}_FlatProc.err"

  cmd="qsub ${qsubopt} -o ${log} -e ${logerr} ./RunThem_wrapper.sh ${file} ${output}"
  #echo $cmd
  eval $cmd

done
