#!/bin/bash

theTime=`date +"%D %T"`
echo $theTime

# /afs/cern.ch/user/a/aconsnd/LaserMeasurements/scripts/createRunList.sh $path $run $events
echo $theTime: runs submitted for first stage analysis >> /afs/cern.ch/user/a/aconsnd/Timing/condorsubmission.log

condor_submit /afs/cern.ch/user/a/aconsnd/LaserMeasurements/scripts/condor_scripts/analyse.sub
