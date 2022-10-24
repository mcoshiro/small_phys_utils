#!/bin/bash

#set CMSSW environment
. /cvmfs/cms.cern.ch/cmsset_default.sh

RUN_KERNEL=$(uname -r | cut -d '-' -f1)

if [ "$RUN_KERNEL" == "3.10.0" ]; then
  export SCRAM_ARCH=slc7_amd64_gcc700
  cd /net/cms29/cms29r0/pico/cc7/CMSSW_10_2_11_patch1/src
elif [ "$RUN_KERNEL" == "2.6.32" ]; then
  cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
fi

. /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -

#set up paths
export SMALL_PHYS_UTILS_DIR=$(dirname $(readlink -e "$BASH_SOURCE"))
export PATH=$SMALL_PHYS_UTILS_DIR/bin::$PATH
export PYTHONPATH=$SMALL_PHYS_UTILS_DIR/libs:$PYTHONPATH
export LD_LIBRARY_PATH=$SMALL_PHYS_UTILS_DIR/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$SMALL_PHYS_UTILS_DIR/inc:$ROOT_INCLUDE_PATH
