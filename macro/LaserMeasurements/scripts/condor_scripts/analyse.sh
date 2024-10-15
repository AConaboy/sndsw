#!/bin/bash

run=$1

source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
echo 'setUp sourced'
source /afs/cern.ch/work/a/aconsnd/HTCondor/config.sh
echo 'config sourced'

printf "Analysing run %06d\n" $run
python3 /afs/cern.ch/user/a/aconsnd/LaserMeasurements/scripts/uproot-analyseRun.py --runNumber $run --make-hists 