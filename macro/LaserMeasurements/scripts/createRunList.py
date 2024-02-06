#!/usr/bin/env python3

def main():
    
    # Clear runlist.txt
    outfile="/afs/cern.ch/user/a/aconsnd/LaserMeasurements/scripts/condor_scripts/runlist.txt"
    with open(outfile, 'w') as f:
        pass  
    
    # Data path
    path='/eos/experiment/sndlhc/raw_data/commissioning/US_tests_LabLausanne_2023/data/'
    
    