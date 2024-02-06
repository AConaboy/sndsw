#!/usr/bin/env python3

from scifi_daq.daq_board.daq_board_multi import DaqBoards
from time import sleep
import logging
import argparse
from datetime import datetime

def main():
  """ Runs the calibration on all boards in the configuration.
  """
  logging.basicConfig(level=logging.INFO)

  parser = argparse.ArgumentParser()
  parser.add_argument('--conf', '-c', type=str, default='conf_example')
  parser.add_argument('--board-ids', '-b', type=int, nargs='+')
  parser.add_argument('--no-tia', action='store_true', help='Skip TIA calibration')
  parser.add_argument('--no-thr', action='store_true', help='Skip THR baseline calibration')
  parser.add_argument('--no-tdc-acq', action='store_true', help='Skip acquisition of TDC calibration data')
  parser.add_argument('--no-tdc-cal', action='store_true', help='Skip calculation of TDS calibration')
  parser.add_argument('--no-qdc-acq', action='store_true', help='Skip acquisition of QDC calibration data')
  parser.add_argument('--no-qdc-cal', action='store_true', help='Skip calculation of QDS calibration')
  parser.add_argument('--no-dcr', action='store_true', help='Skip DCR threshold scan acquisition')
  args = parser.parse_args()


  # the first part of the calibration needs to be performed below the BD voltage
  input("Change HV (set to idle)...")

  d = DaqBoards(args.conf, subset=args.board_ids)
  print("connect")
  start = datetime.now()
  d.connect()
  print('Time elapsed (connect):', datetime.now() - start)

  print("reset")
  start = datetime.now()
  d.fullReset()
  print('Time elapsed (reset):', datetime.now() - start)

  print("initialize")
  start = datetime.now()
  d.initialize()
  
  d.disableLeds(True)
  sleep(1)
  print('Time elapsed (initialize):', datetime.now() - start)

  if not args.no_tia:
    print("TIA baseline calibration")
    start = datetime.now()
    d.tiaCalibration()
    print('Time elapsed (TIA baseline calibration):', datetime.now() - start)

  if not args.no_thr:
    print("THR baseline calibration")
    start = datetime.now()
    d.thrCalibration()
    print('Time elapsed (THR baseline calibration):', datetime.now() - start)

  if not args.no_tdc_acq and not args.no_tdc_cal:
    print("TDC calibration")
    start = datetime.now()
  if not args.no_tdc_acq:
    d.acquireTdcCalibrationData()
  if not args.no_tdc_cal:
    d.calculateTdcCalibration()
  if not args.no_tdc_acq and not args.no_tdc_cal:
    print('Time elapsed (TDC calibration):', datetime.now() - start)


  if not args.no_qdc_acq and not args.no_qdc_cal:
    print("QDC calibration")
  start = datetime.now()
  if not args.no_qdc_acq:
    d.acquireQdcCalibrationData()
  if not args.no_qdc_cal:
    d.calculateQdcCalibration()
  if not args.no_qdc_acq and not args.no_qdc_cal:
    print('Time elapsed (QDC calibration):', datetime.now() - start)

  if not args.no_dcr:
    input("Change HV (set to operation)...")

    print("DCR, thresholds")
    start = datetime.now()
    #set the settings to what you need
    d.acquireThresholdScanData(1, 10, 1) # arguments are: minimum time (per point), minimum number of counts, maximum time
    
    print('Time elapsed (DCR, thresholds):', datetime.now() - start)

if __name__ == '__main__':
  main()