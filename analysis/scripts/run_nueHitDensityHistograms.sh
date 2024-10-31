#!/bin/bash

BDT_flag=$1

if [ -z ${SNDSW_ROOT+x} ];then
	source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
	echo 'setUp sourced'
	source /afs/cern.ch/work/a/aconsnd/HTCondor/config.sh # alienv load sndsw/latest > HTCondor/config.sh
	echo 'config sourced'
fi

if [[ $BDT_flag = true ]]; then
	python $SNDSW_ROOT/analysis/scripts/nueHitDensityHistograms.py --BDTcut
elif [[ $BDT_flag = false ]]; then 
	python $SNDSW_ROOT/analysis/scripts/nueHitDensityHistograms.py
fi
