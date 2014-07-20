#!/usr/bin/env python

import os
import logging
import argparse
from atooms import trajectory

trj_map = {
    'auto': trajectory.Trajectory,
    'xyz' : trajectory.TrajectoryXYZ,
    'rumd': trajectory.TrajectoryRUMD,
    'hoomd': trajectory.TrajectoryHOOMD,
    'lammps': trajectory.TrajectoryLAMMPS,
    'hdf5': trajectory.TrajectoryHDF5,
    }

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--first',   dest='first', type=int, default=0, help='first cfg')
parser.add_argument('-l', '--last',    dest='last', type=int, default=-1, help='last cfg')
parser.add_argument('-s', '--skip',    dest='skip', type=int, default=1, help='interval between cfg')
parser.add_argument('-i', '--fmt-inp', dest='inp', type=str, default='auto', help='input format ')
parser.add_argument('-o', '--fmt-out', dest='out', type=str, default='', help='output format for conversion')
parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
parser.add_argument(nargs='+',         dest='file',type=str, help='input files')
args = parser.parse_args()

if len(args.out) == 0:
    raise ValueError('You must provide an output format')

if not args.out in trj_map:
    raise ValueError('Unknown output format')

if args.out  == 'auto':
    raise ValueError('Cannot use factory for output format')

if args.first == 0 and args.last == -1 and args.skip == 1:
    sl = None
else:
    sl = slice(args.first, args.last, args.skip)

for finp in args.file:
    if not os.path.exists(finp):
        logging.warn('file %s does not exists, skipping it.' % finp)
        continue

    with trj_map[args.inp](finp) as t:
        if sl is None:
            tn = trajectory.NormalizeId(t)
        else:
            tnn = trajectory.NormalizeId(t)
            tn = trajectory.Sliced(t, sl)

        fout = trajectory.convert(tn, trj_map[args.out], args.tag)
        print 'converted %s to %s' % (finp, fout)
    
