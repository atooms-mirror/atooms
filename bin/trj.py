#!/usr/bin/env python

"""Convert trajectory file to a different format."""

import os
import sys
import logging
import argparse
from atooms import trajectory
from atooms.utils import fractional_slice, add_first_last_skip


def print_available_formats():
    print 'Available trajectory formats:'
    fmts = trajectory.Trajectory.formats
    maxlen = max([len(name) for name in fmts])
    for name in sorted(fmts):
        class_name = fmts[name]
        if class_name.__doc__:
            docline = class_name.__doc__.split('\n')[0].rstrip('.')
        else:
            docline = '...no description...'
        fmt = ' - %-' + str(maxlen) + 's : %s'
        print  fmt % (name, docline)

def main(args):
    """Convert trajectory `file_inp` to `file_out`."""
    if args.file_out == '-':
        args.file_out = '/dev/stdout'

    t = trajectory.Trajectory(args.file_inp, fmt=args.inp)

    # If no output format is provided we use the input one
    if args.out is None:
        out_class = t.__class__
    else:
        out_class = args.out

    if args.precision is not None:
        t.precision = args.precision

    if args.flatten_steps:
        t.steps = range(1,len(t)+1)

    # Define slice.
    # We interpret --first N --last N as a request of step N
    if args.last == args.first:
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
        add_interaction_hdf5(fout, args.ff)
    
    if args.file_out != '/dev/stdout':
        print '%s' % fout

    t.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = add_first_last_skip(parser)
    parser.add_argument(      '--fmt-available', dest='fmt_available', action='store_true', help='list available formats')
    parser.add_argument(      '--fmt-fields', dest='fmt', help='format fields')
    parser.add_argument('-I', '--fmt-include', dest='fmt_include', type=str, default='', help='include patterns in format')
    parser.add_argument('-E', '--fmt-exclude', dest='fmt_exclude', type=str, default='', help='exclude patterns from format')
    parser.add_argument('-i', '--fmt-inp', dest='inp', help='input format ')
    parser.add_argument('-o', '--fmt-out', dest='out', help='output format for conversion')
    parser.add_argument('-t', '--tag', dest='tag', type=str, default='', help='tag to add before suffix')
    parser.add_argument('-F', '--ff', dest='ff', type=str, default='', help='force field file')
    parser.add_argument(      '--flatten-steps',dest='flatten_steps', action='store_true', help='use sample index instead of steps')
    parser.add_argument(      '--density', dest='rho', type=float, default=None, help='new density')
    parser.add_argument('-T', '--temperature', dest='temperature', type=float, default=None, help='new temperature')
    parser.add_argument(      '--precision', dest='precision', type=int, default=None, help='write precision')
    parser.add_argument(      '--alphabetic',dest='alphabetic_ids', action='store_true', help='reassign names alphabetically')
    parser.add_argument(nargs=1, dest='file_inp', default='-', help='input file')
    parser.add_argument(nargs='?', dest='file_out', default='-', help='output file')
    args = parser.parse_args()

    if args.fmt_available:
        print_available_formats()
        sys.exit()

    if args.fmt is not None:
        args.fmt = args.fmt.split(',')

    if args.out is not None and not args.out in trajectory.Trajectory.formats:
        print_available_formats()   
        raise ValueError('Unknown output format %s' % args.out)

    args.file_inp = args.file_inp[0]

    main(args)
