#!/usr/bin/env python

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

def add_interaction_hdf5(finp, ff, tag=None):
    """Add interaction to hdf5 file"""
    # TODO: is this a method in hdf5 or not?
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

class Scaling(object):

    def __new__(cls, component, scale, T):
        cls = type('Scaling', (Scaling, component.__class__), component.__dict__)
        return object.__new__(cls)

    def __init__(self, component, factor=None, T=None):
        self.component = component
        self.factor = factor
        self.T = T

    def read_sample(self, *args, **kwargs):
        s = super(Scaling, self).read_sample(*args, **kwargs)

        if self.factor is not None:
            s.cell.side *= self.factor
            for p in s.particle:
                p.position *= self.factor

        if self.T is not None:
            for p in s.particle:
                p.maxwellian(self.T)
            from atooms.system.particle import velocity_cm
            v_cm = velocity_cm(s.particle)
            for p in s.particle:
                p.velocity -= v_cm

        return s

def main(t, args):
    # If no output format is provided we use the input one
    if args.out is None:
        out_class = t.__class__
    else:
        out_class = args.out

    if args.precision is not None:
        t.precision = args.precision

    if args.flatten_steps:
        t.steps = range(1,len(t)+1)
    tn = trajectory.NormalizeId(t)
    t.register_callback(trajectory.decorators.normalize_id, args.alphabetic_ids)

    if args.step:
        step = t.steps.index(args.step)
        ts = trajectory.Sliced(t, slice(step, step+1))
    else:
        # Define slice
        sl = fractional_slice(args.first, args.last, args.skip, len(t))
        # Here we could you a trajectory slice t[sl] but this will load 
        # everything in ram (getitem doesnt provide a generator)
        ts = trajectory.Sliced(t, sl)

    # Change density and
    if args.rho is not None or args.temperature is not None:
        if args.rho is not None:
            rho_0 = t[0].density
            scale = (rho_0 / args.rho)**(1./3)
        tt = Scaling(ts, scale, args.temperature)
    else:
        tt = ts

    fout = trajectory.convert(tt, out_class, tag=args.tag,
                              stdout=args.stdout, fmt=args.fmt,
                              include=args.fmt_include.split(','),
                              exclude=args.fmt_exclude.split(','))

    if args.ff:
        if os.path.exists(args.ff):
            add_interaction_hdf5(fout, args.ff)
        else:
            raise IOError('force field file does not exist')
    
    if not args.stdout:
        print '%s' % fout


parser = argparse.ArgumentParser()
parser = add_first_last_skip(parser)
parser.add_argument(      '--fmt-available', dest='fmt_available', action='store_true', help='list available formats')
parser.add_argument(      '--fmt-fields', dest='fmt', type=str, default=None, help='format fields')
parser.add_argument('-I', '--fmt-include', dest='fmt_include', type=str, default='', help='include patterns in format')
parser.add_argument('-E', '--fmt-exclude', dest='fmt_exclude', type=str, default='', help='exclude patterns from format')
parser.add_argument('-i', '--fmt-inp', dest='inp', type=str, default=None, help='input format ')
parser.add_argument('-o', '--fmt-out', dest='out', type=str, default=None, help='output format for conversion')
parser.add_argument('-S', '--stdout',  dest='stdout', action='store_true', help='dump to stdout')
parser.add_argument('-t', '--tag',     dest='tag', type=str, default='', help='tag to add before suffix')
parser.add_argument('-F', '--ff',      dest='ff', type=str, default='', help='force field file')
parser.add_argument('-g', '--gather',  dest='gather', action='store_true', help='gather files as a super trajectory in <dirname>-conv')
parser.add_argument(      '--flatten-steps',dest='flatten_steps', action='store_true', help='use sample index instead of steps')
parser.add_argument(      '--step',    dest='step', action='store', default=None, type=int, help='')
parser.add_argument(      '--rho',     dest='rho', type=float, default=None, help='new density')
parser.add_argument('-T', '--temperature',dest='temperature', type=float, default=None, help='new temperature')
parser.add_argument(      '--precision',dest='precision', type=int, default=None, help='write precision')
parser.add_argument(      '--alphabetic',dest='alphabetic_ids', action='store_true', help='reassign names alphabetically')
parser.add_argument(nargs='*',         dest='file',type=str, help='input files')
args = parser.parse_args()

if args.fmt_available:
    print_available_formats()
    sys.exit()

if args.fmt is not None:
    args.fmt = args.fmt.split(',')

if len(args.file)==0:
    parser.print_help()

if args.out is not None and not args.out in trajectory.Trajectory.formats:
    print_available_formats()   
    raise ValueError('Unknown output format %s' % args.out)

if args.out  == 'auto':
    raise ValueError('Cannot use factory for output format')

if args.gather:
    with trajectory.SuperTrajectory2(args.file,
                                     trajectory.available_formats[args.inp]) as t:
        main(t, args)
else:
    for finp in args.file:
        if not os.path.exists(finp):
            logging.warn('file %s does not exists, skipping it.' % finp)
            continue
        with trajectory.Trajectory(finp, fmt=args.inp) as t:
            try:
                main(t, args)
            except IOError, e:
                print 'Error: conversion failed for %s (%s)' % (finp, e)
                continue
