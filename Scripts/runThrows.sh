#!/bin/bash
flattreepath="/vols/t2k/users/wparker2/FlatTreesP6T/"
outputpath="/vols/t2k/users/wparker2/Throws2020/"
executable="RunSystBinCorr.exe"
count=0
MAXJOBS=400

run="run2a"
for i in {000..058}
do

    number=$(qstat -u ${USER} | wc -l)
    while [[ ${number} -ge ${MAXJOBS} ]]; do
	number=$(qstat -u ${USER} | wc -l)
	echo "Sleeping to wait for batch jobs to finish..."
	sleep 50
    done

    outdir=${outputpath}/${run}
    infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
    outfile=${outdir}/${run}Throw_${i}.root

    echo "Doing for ${infile}"

    if [ -e ${infile} ]
    then

	if [ -e ${outfile} ]
	then
	    echo "Already done for Run $i"
	else

	    count=1
			
	    scriptDir=${outdir}/scripts
	    logDir=${outdir}/log
	    mkdir -p ${outdir}
	    mkdir -p ${scriptDir}
	    mkdir -p ${logDir}
			
	    ScriptFileName=${scriptDir}/${run}_${i}.sh
	    anaCommand="${executable} -i ${infile} -o ${outfile}"
	    echo "#!/bin/bash" > $ScriptFileName
	    echo "export CMTPATH=/vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3" >> $ScriptFileName
	    echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
	    echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
	    echo "source setup.sh" >> $ScriptFileName
	    echo "cd .." >> $ScriptFileName
	    echo "${anaCommand}" >> $ScriptFileName
	    chmod 744 $ScriptFileName
			
	    stdout="-o ${logDir}/${run}_${i}.log"
	    stderr="-e ${logDir}/${run}_${i}.err"
			
	    qsubmit="qsub"
	    qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
	    command="${ScriptFileName}"
	    #eval ${anaCommand}
#	    eval ${ScriptFileName}
	    eval $qsubmit
			
	fi
    fi
    if (( count == 1 ))
    then
	break
    fi
done

# sleep 300

# run="run2w"
# for i in {000..42}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
#         number=$(qstat -u ${USER} | wc -l)
#         echo "Sleeping to wait for batch jobs to finish..."
#         sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
#         then
#             echo "Already done for Run $i"
#         else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

#         fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run3"
# for i in {000..109}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
#         number=$(qstat -u ${USER} | wc -l)
#         echo "Sleeping to wait for batch jobs to finish..."
#         sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
#         then
#             echo "Already done for Run $i"
#         else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

#         fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run4a"
# for i in {000..128}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run4w"
# for i in {000..128}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run5"
# for i in {000..079}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run6"
# for i in {000..124}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run7"
# for i in {000..118}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run8a"
# for i in {000..128}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run8w"
# for i in {000..096}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done

# sleep 300

# run="run9"
# for i in {000..079}
# do

#     number=$(qstat -u ${USER} | wc -l)
#     while [[ ${number} -ge ${MAXJOBS} ]]; do
# 	number=$(qstat -u ${USER} | wc -l)
# 	echo "Sleeping to wait for batch jobs to finish..."
# 	sleep 50
#     done

#     outdir=${outputpath}/${run}
#     infile=${flattreepath}/${run}/${run}_FlatTree_v2r41_11_25_${i}.root
#     outfile=${outdir}/${run}Throw_${i}.root

#     if [ -e ${infile} ]
#     then

#         if [ -e ${outfile} ]
# 	    then
#             echo "Already done for Run $i"
# 	    else

#             count=1

#             scriptDir=${outdir}/scripts
#             logDir=${outdir}/log
#             mkdir -p ${outdir}
#             mkdir -p ${scriptDir}
#             mkdir -p ${logDir}

#             ScriptFileName=${scriptDir}/${run}_${i}.sh
#             anaCommand="${executable} -i ${infile} -o ${outfile}"
#             echo "#!/bin/bash" > $ScriptFileName
#             echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
#             echo "source /vols/build/t2k/wparker2/MaCh3_ND2019/root/bin/thisroot.sh" >>$ScriptFileName
#             echo "cd /vols/build/t2k/wparker2/MaCh3_ND2019/MaCh3/psyche/psycheSteering/v3r31/cmt" >> $ScriptFileName
#             echo "source setup.sh" >> $ScriptFileName
#             echo "cd .." >> $ScriptFileName
#             echo "${anaCommand}" >> $ScriptFileName
#             chmod 744 $ScriptFileName

#             stdout="-o ${logDir}/${run}_${i}.log"
#             stderr="-e ${logDir}/${run}_${i}.err"

#             qsubmit="qsub"
#             qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
#             command="${ScriptFileName}"
#             #eval ${ScriptFileName}
#             eval $qsubmit

# 	    fi
#     fi
# #    if (( count == 1 ))
# #    then
# #        break
# #    fi
# done
