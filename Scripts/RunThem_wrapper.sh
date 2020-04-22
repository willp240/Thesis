#!/bin/bash

if [[ $# -ne 2 ]]; then
  echo "$0 needs two arguments, you gave me $#"
  echo "$0 InputFile -o Outputfile"
  exit
fi

SLEEPMAX=60
PYSLEEP=$(python -c "from random import random; print int(random()*${SLEEPMAX});")
echo "Sleeping ${PYSLEEP}"
sleep ${PYSLEEP}

./RunSystBinCorr.exe -i $1 -o $2
