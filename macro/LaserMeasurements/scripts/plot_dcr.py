#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import os.path as op

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('conf', type=str, help='path to the BOARD configuration folder (e.g. conf_example/board_x)')
  args = parser.parse_args()
  
  with h5py.File(op.join(args.conf, 'pe_cal/data.hdf5'), 'r') as f:
    for i in [2, 3]:
      fig1, ax1 = plt.subplots()
      fig2, ax2 = plt.subplots()
      fig3, ax3 = plt.subplots()

      ax1.plot(np.array(f[f'T1/tofpet_{i}/rate']).T)
      ax1.set_title(f'tofpet_{i}, T1')
      ax1.set_xlabel('Threshold T1')
      ax1.set_ylabel('Rate [Hz]')
      ax1.grid(alpha=0.3)
      ax1.set_yscale('log')

      ax2.plot(np.array(f[f'T2/tofpet_{i}/rate']).T)
      ax2.set_title(f'tofpet_{i}, T2')
      ax2.set_xlabel('Threshold T2')
      ax2.set_ylabel('Rate [Hz]')
      ax2.grid(alpha=0.3)
      ax2.set_yscale('log')

      ax3.plot(np.array(f[f'E/tofpet_{i}/rate']).T)
      ax3.set_title(f'tofpet_{i}, E')
      ax3.set_xlabel('Threshold E')
      ax3.set_ylabel('Rate [Hz]')
      ax3.grid(alpha=0.3)
      ax3.set_yscale('log')

      plt.show()

if __name__ == '__main__':
  main()