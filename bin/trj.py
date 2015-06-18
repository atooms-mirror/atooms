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

def add_interaction_hdf5(finp, ff, tag=None):
    """Add interaction to hdf5 file"""

    import h5py

    pid = os.getpid()
    f_ref = '/tmp/cnv_%s.h5' % pid
    # TODO: we can cache a ref file if ff is the same
    os.system('system.x -n 2 -f %s %s > /dev/stdout' % (ff, f_ref))
    ref = h5py.File(f_ref, 'r')

    if tag:
        d = os.path.dirname(finp) + '_' + opts.tag
        if not os.path.exists(d):
            os.makedirs(d)
        fout = d + '/' + os.path.basename(finp)
    else:
        fout = finp + '.bak'
        
    # Add interaction
    os.system('/bin/cp %s %s' % (finp, fout))
    h5 = h5py.File(fout , 'r+')
    # Make sure interaction does not exist
    try:
        del h5['initialstate/interaction']
    except:
        pass
    h5.copy(ref['initialstate/interaction'], 'initialstate/interaction')
    h5.close()

    if tag is None:
        os.remove(finp)
        os.rename(fout, finp)

    # Final cleanup
    ref.close()
    import glob
    for f in glob.glob(f_ref + '*'):
        os.remove(f)


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--percentage',dest='percent', action='store_true', help='first and last in percentages (if set)')
parser.add_argument('-f', '--first',   dest='first', type=int, default=0, help='first cfg')
parser.add_argument('-l', '--last',    dest='last', type=int, default=-1, help='last cfg')
parser.add_argument('-s', '--skip',    dest='skip', type=int, default=1, help='interval between cfg')
parser.add_argument('-i', '--fmt-inp', dest='inp', type=str, default='auto', help='input format ')
parser.add_argument('-o', '--fmt-out', dest='out', type=str, default='', help='output format for conversion')
parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
parser.add_argument('-F', '--ff',      dest='ff', type=str, default='', help='force field file')
parser.add_argument(nargs='+',         dest='file',type=str, help='input files')
args = parser.parse_args()

if len(args.out) == 0:
    raise ValueError('You must provide an output format')

if not args.out in trj_map:
    raise ValueError('Unknown output format')

if args.out  == 'auto':
    raise ValueError('Cannot use factory for output format')

for finp in args.file:
    if not os.path.exists(finp):
        logging.warn('file %s does not exists, skipping it.' % finp)
        continue

    with trj_map[args.inp](finp) as t:
        tn = trajectory.NormalizeId(t)

        # Define slice
        if args.percent:
            # If we ask fractional first / last sample and they are different from default
            # get the actual first and/or last sample to convert
            if args.first > 0:
                args.first = int(args.first / 100. * len(tn))
            if args.last != -1:
                args.last = int(args.last / 100. * len(tn))

        sl = slice(args.first, args.last, args.skip)
        t = trajectory.Sliced(t, sl)
        fout = trajectory.convert(t, trj_map[args.out], args.tag)
        print 'converted %s to %s from step %d to %d' % (finp, fout, args.first, args.last)
    
        if args.ff:
            if os.path.exists(args.ff):
                add_interaction_hdf5(fout, args.ff)
            else:
                raise OSError('force field file does not exist')
