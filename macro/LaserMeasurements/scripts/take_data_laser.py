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

    print('Time elapsed (initializing):', datetime.now() - start)

    nTriggers = args.runtime * 1000

    # for ph in range(0, 224, 1):
    for ph in [0]:
      d.daqBoards.setTrigger(sources=['sequencer'], seqCounter=nTriggers, seqPeriod=40000)
      d.daqBoards.setInjection(enable=True, phase=ph, duration=1, clk='clk40')
      d.daqBoards.setDebugBoardOutput(sma0='test_pulse', sma1='trigger')

      print("starting")
      start = datetime.now()
      d.startDaq(args.run if args.run >= 0 else None, saveConfiguration=True)
      d.daqBoards.startTriggerSequencer(False)
      print('Time elapsed (starting):', datetime.now() - start)

      start = datetime.now()

      runNumber = d.daqServer.getStatus()['run_number']
      print(f'Run {runNumber}, {conf_dir}, phase {ph}')
      with open(f'logs/log_{runNumber}.txt', 'w') as f:
        with tqdm(total=nTriggers, unit=' triggers') as t:
          while True:
            status = d.daqBoards.getTrigger()[26]
            t.update(nTriggers - status['seq_counter_remaining'] - t.n)
            f.write(str(datetime.now()) + '\n')
            f.write(str(d.daqServer.getStatus()) + '\n')
            f.write(str(d.daqBoards.getStatus()) + '\n')
            f.flush()
            if status['seq_counter_remaining'] == 0:
              break
            sleep(1)

      # stops everything
      print(d.stopDaq())
      print('Time elapsed (DAQ):', datetime.now() - start)
      return args
      
if __name__ == '__main__':
  args=main()
  
  # automatically produce hists after data taking
  # options.make_hists=True
  # ld=LaserData(args)
  # ld.EventLoop()
  # ld.WriteOutHistograms()      


