#!/bin/bash

# geo_file="/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root"
geo_file="/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/geofile_full.PG-TGeant4-0.root"

input_file=$1
outputfile=$2
mode=$3
cut_set=$4

if [ -z ${SNDSW_ROOT+x} ];then
	source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
	echo 'setUp sourced'
	source /afs/cern.ch/work/a/aconsnd/HTCondor/config.sh # alienv load sndsw/latest > HTCondor/config.sh
	echo 'config sourced'
fi

mode_list=(stage1cuts novetocuts FVsideband allowWalls2and5 stage1cutsVetoFirst nueFilter)

if [[ ! ${mode_list[@]} =~ $MODE ]]; then
    echo Mode $mode not available. It must be one of "${mode_list[*]}"
    echo Exitting.
    exit
fi

echo cut set $cut_set for mode $mode

if [ "$mode" == "nueFilter" ]; then 

	# Run first stage filter
	neutrinoFilterGoldenSample $input_file $outputfile $cut_set
fi