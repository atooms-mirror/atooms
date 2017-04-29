#!/usr/bin/env python

"""Convert trajectory file to a different format."""

import os
import sys
import logging
import argparse
import random
from atooms import trajectory
from atooms.utils import fractional_slice, add_first_last_skip


def available_formats():
    txt = 'available trajectory formats:\n'
    fmts = trajectory.Trajectory.formats
    maxlen = max([len(name) for name in fmts])
    for name in sorted(fmts):
        class_name = fmts[name]
        if class_name.__doc__:
            docline = class_name.__doc__.split('\n')[0].rstrip('.')
        else:
            docline = '...no description...'
        fmt = '  %-' + str(maxlen) + 's : %s\n'
        txt += fmt % (name, docline)
    return txt

def info(trajectory):
    from atooms.system.particle import species, composition
    txt = ''
    txt += 'path                 = %s\n' % trajectory.filename
    txt += 'format               = %s\n' % trajectory.__class__
    txt += 'frames               = %s\n' % len(trajectory)
    txt += 'megabytes            = %s\n' % (os.path.getsize(trajectory.filename) / 1e6)
    txt += 'particles            = %s\n' % len(trajectory[0].particle)
    txt += 'species              = %s\n' % len(species(trajectory[0].particle))
    txt += 'composition          = %s\n' % list(composition(trajectory[0].particle))
    txt += 'density              = %s\n' % round(trajectory[0].density, 10)
    txt += 'cell side            = %s\n' % trajectory[0].cell.side
    txt += 'cell volume          = %s\n' % trajectory[0].cell.volume
    if len(trajectory)>1:
        txt += 'steps                = %s\n' % trajectory.steps[-1]
        txt += 'duration             = %s\n' % trajectory.times[-1]
        txt += 'timestep             = %s\n' % trajectory.timestep
        txt += 'block size           = %s\n' % trajectory.block_size
        if trajectory.block_size == 1:
            txt += 'steps between frames = %s\n' % (trajectory.steps[1]-trajectory.steps[0])
            txt += 'time between frames  = %s\n' % (trajectory.times[1]-trajectory.times[0])
        else:
            txt += 'block steps          = %s\n' % trajectory.steps[trajectory.block_size]
            txt += 'block                = %s\n' % ([trajectory.steps[i] for i in range(trajectory.block_size+1)])
        txt += 'grandcanonical       = %s' % trajectory.grandcanonical
    print txt

def main(args):
    """Convert trajectory `file_inp` to `file_out`."""
    if args.file_out == '-':
        args.file_out = '/dev/stdout'

    if args.folder:
        t = trajectory.folder.Foldered(args.file_inp, cls=args.inp)
    else:
        t = trajectory.Trajectory(args.file_inp, fmt=args.inp)

    if args.info:
        info(t)
        return

    # If no output format is provided we use the input one
    if args.out is None:
        out_class = t.__class__
    else:
        out_class = args.out

    if args.precision is not None:
        t.precision = args.precision

    if args.flatten_steps:
        t.steps = range(1,len(t)+1)

    # Reset random number generator
    if args.seed:
        random.seed(args.seed)

    # Trick to allow some trajectory formats to set the box side.
    # This way the cell is defined as we read the sample (callbacks
    # will not do that).
    if args.side is not None:
        t._side = args.side

    # Define slice.
    # We interpret --first N --last N as a request of step N
    if args.last == args.first and args.last is not None:
        args.last += 1
    sl = fractional_slice(args.first, args.last, args.skip, len(t))
    # Here we could you a trajectory slice t[sl] but this will load
    # everything in ram (getitem doesnt provide a generator). This
    # will be fixed with python 3.
    ts = trajectory.Sliced(t, sl)

    # Change density and temperature
    if args.rho is not None:
        ts.register_callback(trajectory.decorators.set_density, args.rho)
    if args.temperature is not None:
        ts.register_callback(trajectory.decorators.set_temperature, args.temperature)

    # We always normalize species id's using fortran convention
    ts.register_callback(trajectory.decorators.normalize_id, args.alphabetic_ids)

    # Trajectory conversion
    fout = trajectory.convert(ts, out_class, args.file_out,
                              tag=args.tag, fmt=args.fmt,
                              include=args.fmt_include.split(','),
                              exclude=args.fmt_exclude.split(','))

    if args.ff:
        from atooms.trajectory.hdf5 import add_interaction_hdf5
        add_interaction_hdf5(fout, args.ff)
    
    if args.file_out != '/dev/stdout':
        print '%s' % fout

    t.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(epilog=available_formats(), 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser = add_first_last_skip(parser)
    parser.add_argument(      '--fmt-fields', dest='fmt', help='format fields')
    parser.add_argument('-I', '--fmt-include', dest='fmt_include', type=str, default='', help='include patterns in format')
    parser.add_argument('-E', '--fmt-exclude', dest='fmt_exclude', type=str, default='', help='exclude patterns from format')
    parser.add_argument('-i', '--fmt-inp', dest='inp', help='input format ')
    parser.add_argument('-o', '--fmt-out', dest='out', help='output format for conversion')
    parser.add_argument(      '--folder', dest='folder', action='store_true', help='force folder-based layout')
    parser.add_argument('-t', '--tag', dest='tag', type=str, default='', help='tag to add before suffix')
    parser.add_argument('-F', '--ff', dest='ff', type=str, default='', help='force field file')
    parser.add_argument(      '--flatten-steps',dest='flatten_steps', action='store_true', help='use sample index instead of steps')
    parser.add_argument(      '--side', dest='side', type=float, default=None, help='set cell side')
    parser.add_argument(      '--density', dest='rho', type=float, default=None, help='new density')
    parser.add_argument('-T', '--temperature', dest='temperature', type=float, default=None, help='new temperature')
    parser.add_argument(      '--precision', dest='precision', type=int, default=None, help='write precision')
    parser.add_argument(      '--alphabetic',dest='alphabetic_ids', action='store_true', help='reassign names alphabetically')
    parser.add_argument(      '--info', dest='info', action='store_true', help='print info')
    parser.add_argument(      '--seed', dest='seed', type=int, help='set seed of random number generator')
    parser.add_argument(nargs=1, dest='file_inp', default='-', help='input file')
    parser.add_argument(nargs='?', dest='file_out', default='-', help='output file')
    args = parser.parse_args()

    if args.fmt is not None:
        args.fmt = args.fmt.split(',')

    if args.out is not None and not args.out in trajectory.Trajectory.formats:
        available_formats()   
        raise ValueError('Unknown output format %s' % args.out)

    args.file_inp = args.file_inp[0]

    main(args)
