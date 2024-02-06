#!/usr/bin/env python3

from scifi_daq.daq import Daq
from time import sleep
import logging
import argparse
from datetime import datetime
from tqdm import tqdm

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('run', type=int, help='Set to -1 to automatically select the run number')
  parser.add_argument('--conf', '-c', default='conf_example')
  parser.add_argument('--board-ids', '-b', type=int, nargs='+')
  parser.add_argument('--runtime', '-t', type=int, default=None)
  parser.add_argument('--events', '-e', type=int, default=20)
  args = parser.parse_args()
  conf_dir = args.conf

  countEvents = args.runtime is None
  totalCounts = args.events if countEvents else args.runtime

  logging.basicConfig(level=logging.ERROR)

  d = Daq(conf_dir, subset=args.board_ids)

  print("initialize")
  start = datetime.now()
  d.connect()
  d.initialize(args.run if args.run >= 0 else None, syncWithOrbit=False)
  d.daqBoards.disableLeds(True)

  sleep(1)

  print('Time elapsed (initializing):', datetime.now() - start)

  nTriggers = args.runtime * 1000

  d.daqBoards.setTrigger(sources=['sequencer'], seqCounter=nTriggers, seqPeriod=40000)
  d.daqBoards.setInjection(enable=True, phase=ph, duration=1, clk='clk40')
  d.daqBoards.setDebugBoardOutput(sma0='test_pulse', sma1='trigger')

  # sets up the DAQ boards to accept triggers from the TTC system
  # d.daqBoards.setTrigger(sources=['ttc'])

  #sets up the TTC system to send triggers when a command is received (triggers are used to check synchronization)
  # d.vmeClient.setL1aTrigger(l1a_input='vme')

  # print("starting")
  # start = datetime.now()
  # d.startDaq(args.run if args.run >= 0 else None, saveConfiguration=True)
  # # d.daqBoards.startTriggerSequencer(False)
  # print('Time elapsed (starting):', datetime.now() - start)

  # start = datetime.now()

  # runNumber = d.daqServer.getStatus()['run_number']
  # print(f'Run {runNumber}')
  # TODO add scan in HV
  # TODO add scan in parameters

  with open('logs/log_daq.txt', 'a') as glog:
    try:
      lastIteration = 0
      while lastIteration < totalCounts - 1:
        print("starting")
        start = datetime.now()
        d.startDaq(runNumber=(args.run if args.run >= 0 else None), saveConfiguration=True)
        print('Time elapsed (starting):', datetime.now() - start)

        start = datetime.now()
        runNumber = d.daqServer.getStatus()['run_number']
        print(f'Starting run: {runNumber}')
        glog.write(f'{datetime.now()} starting run {runNumber}\n')
      
        with open(f'logs/log_{runNumber}.txt', 'w') as f:
          # for i in tqdm(range(totalIterations - lastIteration)):
          with tqdm(total=totalCounts - lastIteration, unit=' events' if countEvents else ' s') as t:
            while True:
              d.vmeClient.generateSoftL1A()
              f.write(str(datetime.now()) + '\n')
              sstatus = d.daqServer.getStatus()
              f.write(str(sstatus) + '\n')
              f.write(str(d.daqBoards.getStatus()) + '\n')

              if countEvents:
                lastIteration = sstatus['daq_monitor']['event_processors_monitor']['processed_events']
              else:
                lastIteration += 1
              
              t.update(lastIteration - t.n)
              
              if lastIteration >= ((args.events) if countEvents else (args.runtime)):
                break
              
              if not sstatus['daq_monitor']['event_builder_trigless_monitor']['sync']:
                print('Desync!')
                glog.write(f'{datetime.now()} desync found at iteration {lastIteration}\n')
                break
              sleep(1)
          
          f.write(str(datetime.now()) + '\n')
          finalstatus = d.stopDaq()
          f.write(str(finalstatus) + '\n')
          f.write(str(d.daqBoards.getStatus()) + '\n')
          print('Time elapsed (DAQ):', datetime.now() - start)
        glog.write(f'{datetime.now()} run {runNumber} finished, events written: {finalstatus["daq_monitor"]["event_writer_monitor"]["written_events"]}\n')
    except KeyboardInterrupt:
      print('Stopping DAQ...')
      finalstatus = d.stopDaq()
      glog.write(f'{datetime.now()} run {runNumber} stopped, events written: {finalstatus["daq_monitor"]["event_writer_monitor"]["written_events"]}\n')
      print(finalstatus)
      print('Time elapsed (DAQ):', datetime.now() - start)


if __name__ == '__main__':
  main()
