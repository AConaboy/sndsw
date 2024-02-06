#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scifi_daq.daq import Daq
from time import sleep
import logging
import argparse
from datetime import date, datetime
from tqdm import tqdm 
# from analyseRun import LaserData

def main():
  # The run number is taken as a command line argument
  parser = argparse.ArgumentParser()
  parser.add_argument('run', type=int)
  parser.add_argument('--runtime', '-t', type=int, default=60)
  parser.add_argument('--mode', '-m', default='qdc')

  # For analyseRun
  parser.add_argument('--runNumber', '-r', dest='runNumber',type=int, default=-1, required=False)
  parser.add_argument('--PCB', dest='PCB',type=str, default='US', required=False)
  parser.add_argument('--make-hists', action='store_true', required=False)

  args = parser.parse_args()

  logging.basicConfig(level=logging.ERROR)

  for conf_dir in ['conf_example']:
#   for conf_dir in ['conf_darkRoom_qdc_off_2302', 'conf_darkRoom_qdc_3ns_2302', 'conf_darkRoom_qdc_6ns_2302', 'conf_darkRoom_tot_off_2302', 'conf_darkRoom_tot_3ns_2302', 'conf_darkRoom_tot_6ns_2302']:
  # for conf_dir in ['conf_darkRoom_qdc_off_2302', 'conf_darkRoom_qdc_3ns_2302', 'conf_darkRoom_qdc_6ns_2302']:

    d = Daq(conf_dir, subset=None)

    print("initialize")
    start = datetime.now()
    d.connect()
    d.initialize(args.run if args.run >= 0 else None)
    d.daqBoards.disableLeds(True)
    sleep(1)

    start = datetime.now()
    d.startDaq(args.run if args.run >= 0 else None, saveConfiguration=True)
    print('Time elapsed (starting):', datetime.now() - start)
    runNumber = d.daqServer.getStatus()['run_number']

    # for ph in range(0, 224, 1):
    phases=[10]
    # phases=range(0,224,1)
    for ph in phases:

      # Same total triggers as internal triggering mode.
      nTriggers = int(args.runtime * 2000 / len(phases))
      
      d.daqBoards.setTrigger(sources=['sequencer'], seqCounter=nTriggers, seqPeriod=40000)
      d.daqBoards.setInjection(enable=True, phase=ph, duration=1, clk='clk40')
      d.daqBoards.setDebugBoardOutput(sma0='test_pulse', sma1='trigger')    
      
      d.daqBoards.startTriggerSequencer(True)
      
      print(f'Run {runNumber}, {conf_dir}, phase {ph}')
      with open(f'logs/log_{runNumber}.txt', 'w') as f:
        
        f.write(f'Trigger phase: {ph} \n')
        f.flush()

    # stops everything
    print(d.stopDaq())
    print(f'All trigger phases cycled.')
    print('Time elapsed (DAQ):', datetime.now() - start)
    return args
      
if __name__ == '__main__':
  args=main()
  