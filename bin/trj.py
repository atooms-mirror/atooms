#!/usr/bin/env python

import os
import logging
import argparse
from atooms import trajectory
from atooms.utils import fractional_slice, add_first_last_skip

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
    os.system('system.x -n 2 -f %s %s 1>/dev/null 2>/dev/null' % (ff, f_ref))
    ref = h5py.File(f_ref, 'r')

    # TODO: can this tag de bropped?
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
parser = add_first_last_skip(parser)
parser.add_argument('-V', '--velocity',dest='vel', action='store_true', help='dump velocities')
parser.add_argument('-i', '--fmt-inp', dest='inp', type=str, default='auto', help='input format ')
parser.add_argument('-o', '--fmt-out', dest='out', type=str, default='', help='output format for conversion')
parser.add_argument('-S', '--stdout',  dest='stdout', action='store_true', help='dump to stdout')
parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
parser.add_argument('-F', '--ff',      dest='ff', type=str, default='', help='force field file')
parser.add_argument(      '--flatten-steps',dest='flatten_steps', action='store_true', help='use sample index instead of steps')
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
        if args.flatten_steps:
            t.steps = range(1,len(t)+1)
        tn = trajectory.NormalizeId(t)

        # Define slice
        sl = fractional_slice(args.first, args.last, args.skip, len(tn))
        # Here we could you a trajectory slice t[sl] but this will load 
        # everything in ram (getitem doesnt provide a generator)
        ts = trajectory.Sliced(tn, sl)
        try:
            if args.vel:
                fout = trajectory.convert(ts, trj_map[args.out], tag=args.tag, stdout=args.stdout, include=['vx','vy','vz'])
            else:
                fout = trajectory.convert(ts, trj_map[args.out], tag=args.tag, stdout=args.stdout)
        except IOError, e:
            print 'Error: conversion failed for %s (%s)' % (finp, e)
            continue

        if args.ff:
            if os.path.exists(args.ff):
                add_interaction_hdf5(fout, args.ff)
            else:
                raise IOError('force field file does not exist')

        print '%s' % fout
