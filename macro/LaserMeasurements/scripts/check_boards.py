#!/usr/bin/env python3

from scifi_daq.daq import Daq
import argparse
import socket

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--conf', '-c', default='conf_example')
  parser.add_argument('--names', action='store_true', help='print names and not IDs')
  parser.add_argument('--board-ids', '-b', type=int, nargs='+')
  args = parser.parse_args()

  d = Daq(args.conf, subset=args.board_ids)

  unreachable_boards = []
  unreachable_boards_names = []

  for n, b in d.daqBoards._boards.items():
    try:
      b.connect(timeout=1)
    except (TimeoutError, socket.timeout):
      unreachable_boards.append(n)
      unreachable_boards_names.append(b.name)
  
  if args.names:
    print(' '.join(map(str, unreachable_boards_names)))
  else:
    print(','.join(map(str, unreachable_boards)))


if __name__ == '__main__':
  main()