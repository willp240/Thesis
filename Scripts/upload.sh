#!/bin/bash

# Loop over the directories, upload
icd /QMULZone1/home/asg/asg2019oa/ND280/Highlandv2r39_TestFlatTrees_18Nov2019

here=$(pwd -P)
there=$(ipwd)
for i in $(ls -d run*); do
  cd ${here}
  icd ${there}
  echo "Uploading files in ${i}/root to ${i} on irods..."
  imkdir -p ${i}
  iput -b -f -v -a -r ${here}/${i}/root ${i}
done
