#!/bin/bash

# Clear runlist.txt
outfile="/afs/cern.ch/user/a/aconsnd/LaserMeasurements/scripts/condor_scripts/runlist.txt"
> "$outfile"

# Data path
path=/eos/experiment/sndlhc/raw_data/commissioning/US_tests_LabLausanne_2023/data/

# Count number of runs
count=$(find "$path" -type d -name "run_*" | wc -l)

# Number of runs to be analysed per job
RunsPerJobMax=10

# printf "%srun_%06d.root\n" "$path" "$count" >> "$outfile"
echo $count $(($count-$RunsPerJobMax)) >> "$outfile"

while ((count >= $RunsPerJobMax)); do
    result=$((count - $RunsPerJobMax))
    # printf "%srun_%06d.root\n" "$path" "$result" >> "$outfile"
    echo $result $(($result - $RunsPerJobMax)) >> "$outfile"
    count=$result
done
