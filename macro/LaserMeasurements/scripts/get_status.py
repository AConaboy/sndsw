#!/usr/bin/env python3

import argparse
from scifi_daq.daq_board.daq_board_multi import DaqBoards
import logging

def main():
  logging.basicConfig(level=logging.INFO)

  parser = argparse.ArgumentParser()
  parser.add_argument('--conf', '-c', type=str, default='conf_example')
  parser.add_argument('--board-ids', '-b', type=int, nargs='+')
  args = parser.parse_args()

  d = DaqBoards(args.conf, subset=args.board_ids)
  d.connect()
  st = d.getStatus()
  for bid, status in st.items():
    print(bid)
    print(status)


if __name__ == '__main__':
  main()