export T2KREWEIGHT=/vols/build/t2k/wparker2/T2KReWeight
export LD_LIBRARY_PATH=$T2KREWEIGHT/lib:$LD_LIBRARY_PATH
export PATH=$T2KREWEIGHT/bin:$PATH

#Root install
source /vols/build/t2k/wparker2/ROOT/v5r34p34n00/Linux-x86_64/bin/thisroot.sh
export PATH=$T2KREWEIGHT/bin:$PATH:$T2KREWEIGHT/app:$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$T2KREWEIGHT/lib:$LD_LIBRARY_PATH;

# Set up analysis reader
here=$(pwd -P)
export CMTPATH=/vols/build/t2k/wparker2/
cd ${CMTPATH}/nd280Highland2/v2r29/cmt
source setup.sh
cd ${here}

# NEUT and CERNLIB
export CERN=/vols/build/t2k/cvw09/CERNLIB
export CERN_LEVEL=2005
export CERN_ROOT=$CERN/$CERN_LEVEL
export CERNLIB=$CERN_ROOT/lib
export LD_LIBRARY_PATH=$CERNLIB:$LD_LIBRARY_PATH
export PATH=$CERN_ROOT/bin:$PATH
export NEUT_ROOT=/vols/build/t2k/wparker2/neut_5.3.3_maqefix
export PATH=$NEUT_ROOT/src/neutsmpl/bin:$PATH
export LD_LIBRARY_PATH=$NEUT_ROOT/src/reweight:$LD_LIBRARY_PATH

# NIWG
export NIWG=/vols/build/t2k/wparker2/NIWGReWeight
export LD_LIBRARY_PATH=$NIWG:$LD_LIBRARY_PATH
export NIWGREWEIGHT_INPUTS=$NIWG/inputs

#JNU 
export JNUBEAM=$CMTPATH/JReWeight/
export LD_LIBRARY_PATH=${JNUBEAM}:$LD_LIBRARY_PATH;
export JREWEIGHT_INPUTS=${JNUBEAM}/inputs

#psyche
export PSYCHELIBS=/vols/build/t2k/wparker2/T2KReWeight/psychelibs/
export PSYCHEINCLUDES=/vols/build/t2k/wparker2/T2KReWeight/psycheinc/
export LD_LIBRARY_PATH=$PSYCHELIBS:$LD_LIBRARY_PATH
