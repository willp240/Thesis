#!/bin/bash

export CMTPATH=/vols/build/t2k/wparker2/

# run="2a"
# for i in {000..184}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_nu.root"
    
#     if [ -e "${infile}" ]
#     then
    
# 	outfile="/vols/t2k/users/wparker2/Throws/Run${run}_TEST/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"
		    
# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}_TEST"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
# 	mkdir -p $scriptDir
# 	mkdir -p $logDir
	
# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
# 	echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
# 	echo "${anaCommand}" >> $ScriptFileName
# 	chmod 744 $ScriptFileName
	
# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"
	
# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
# 	eval $qsubmit
	
#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="2w"
# for i in {000..240}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="3"
# for i in {000..615}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="3b"
# for i in {000..89}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="3c"
# for i in {000..526}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="4a"
# for i in {000..699}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="4w"
# for i in {000..452}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="5"
# for i in {000..419}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="6b"
# for i in {000..715}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="6c"
# for i in {000..266}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="6d"
# for i in {000..347}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="6e"
# for i in {000..433}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6B_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="7"
# for i in {000..674}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6L_antinu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="8a"
# for i in {000..725}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6L_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

# run="8w"
# for i in {000..528}
# do
#     infile="/vols/t2k/users/wparker2/FlatTrees/Run${run}/FlatTree_Run${run}_mc.${i}_prod6L_nu.root"

#     if [ -e "${infile}" ]
#     then

#         outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
# 	echo "Doing for ${infile}"

# 	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
# 	scriptDir=${ScratchDir}/scripts
# 	logDir=${ScratchDir}/log
#         mkdir -p $scriptDir
# 	mkdir -p $logDir

# 	ScriptFileName=${scriptDir}/run${run}_${i}.sh
# 	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
# 	echo "#!/bin/bash" > $ScriptFileName
# 	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
# 	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
#         echo "source setup.sh" >> $ScriptFileName
# 	echo "cd .." >> $ScriptFileName
#         echo "${anaCommand}" >> $ScriptFileName
#         chmod 744 $ScriptFileName

# 	stdout="-o ${logDir}/Run_${i}.log"
# 	stderr="-e ${logDir}/Run_${i}.err"

# 	qsubmit="qsub"
# 	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
# 	command="${ScriptFileName}"
# 	#eval ${ScriptFileName}
#         eval $qsubmit

#     else
# 	echo "${infile} doesn't exist"
#     fi
# done

run="Sand_FHC"
for i in {000..476}
do
    infile="/vols/t2k/users/wparker2/FlatTrees/${run}/FlatTree_sand_mc.${i}_prod6B_nu.root"

    if [ -e "${infile}" ]
    then

        outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
	echo "Doing for ${infile}"

	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
	scriptDir=${ScratchDir}/scripts
	logDir=${ScratchDir}/log
        mkdir -p $scriptDir
	mkdir -p $logDir

	ScriptFileName=${scriptDir}/run${run}_${i}.sh
	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
	echo "#!/bin/bash" > $ScriptFileName
	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
        echo "source setup.sh" >> $ScriptFileName
	echo "cd .." >> $ScriptFileName
        echo "${anaCommand}" >> $ScriptFileName
        chmod 744 $ScriptFileName

	stdout="-o ${logDir}/Run_${i}.log"
	stderr="-e ${logDir}/Run_${i}.err"

	qsubmit="qsub"
	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
	command="${ScriptFileName}"
	#eval ${ScriptFileName}
        eval $qsubmit

    else
	echo "${infile} doesn't exist"
    fi
done

run="Sand_RHC"
for i in {000..471}
do
    infile="/vols/t2k/users/wparker2/FlatTrees/${run}/FlatTree_sand_mc.${i}_prod6B_antinu.root"

    if [ -e "${infile}" ]
    then

        outfile="/vols/t2k/users/wparker2/Throws/Run${run}/Run${run}Throw_${i}.root"
	echo "Doing for ${infile}"

	ScratchDir="/vols/t2k/users/wparker2/Throws/Run${run}"
	scriptDir=${ScratchDir}/scripts
	logDir=${ScratchDir}/log
        mkdir -p $scriptDir
	mkdir -p $logDir

	ScriptFileName=${scriptDir}/run${run}_${i}.sh
	anaCommand="RunSystBinCorr.exe -i ${infile} -o ${outfile}"
	echo "#!/bin/bash" > $ScriptFileName
	echo "export CMTPATH=/vols/build/t2k/wparker2" >> $ScriptFileName
	echo "cd ${CMTPATH}/psyche/psycheSteering/v3r25/cmt" >> $ScriptFileName
        echo "source setup.sh" >> $ScriptFileName
	echo "cd .." >> $ScriptFileName
        echo "${anaCommand}" >> $ScriptFileName
        chmod 744 $ScriptFileName

	stdout="-o ${logDir}/Run_${i}.log"
	stderr="-e ${logDir}/Run_${i}.err"

	qsubmit="qsub"
	qsubmit="${qsubmit} -l h_vmem=16G ${stdout} ${stderr} ${ScriptFileName}"
	command="${ScriptFileName}"
	#eval ${ScriptFileName}
        eval $qsubmit

    else
	echo "${infile} doesn't exist"
    fi
done
