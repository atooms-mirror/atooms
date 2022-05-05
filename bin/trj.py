#!/usr/bin/env python

"""Convert trajectory file to a different format."""

from __future__ import print_function
import os
import sys
import logging
import argparse
import random
from atooms import trajectory
from atooms.core.utils import fractional_slice, add_first_last_skip, setup_logging
from atooms.trajectory.utils import check_block_size, info, formats


def main_info(args):
    """Print info on trajectories."""
    for file_inp in args.file_inp:
        if args.folder:
            th = trajectory.folder.Foldered(file_inp, cls=args.inp)
        else:
            th = trajectory.Trajectory(file_inp, fmt=args.inp)

        if args.fields:
            print(info(th, args.fields))
        else:
            print(info(th))

def main(args):
    """Convert trajectory `file_inp` to `file_out`."""
    args.file_inp = args.file_inp[0]

    if args.fields is not None:
        args.fields = args.fields.split(',')

    if args.out is not None and not args.out in trajectory.Trajectory.formats:
        formats()
        raise ValueError('Unknown output format %s' % args.out)

    if args.file_out == '-':
        args.file_out = '/dev/stdout'

    if args.verbose:
        setup_logging('atooms', 20)
        import atooms.core.progress
        if args.file_out != '/dev/stdout':
            atooms.core.progress.active = True

    if args.folder:
        t = trajectory.folder.Foldered(args.file_inp, cls=args.inp)
    else:
        t = trajectory.Trajectory(args.file_inp, fmt=args.inp)

    # If no output format is provided we use the input one
    if args.out is None:
        out_class = t.__class__
    else:
        out_class = args.out

    if args.precision is not None:
        t.precision = args.precision

    if args.flatten_steps:
        t.steps = range(1, len(t)+1)

    # Reset random number generator
    if args.seed:
        random.seed(args.seed)

    # Trick to allow some trajectory formats to set the box side.
    # This way the cell is defined as we read the sample (callbacks
    # will not do that).
    if args.side is not None:
        def fix_cell(system, side):
            from atooms.system import Cell
            system.cell = Cell(side)
            return system
        L = args.side.split(',')
        if len(L) == 1:
            t.add_callback(fix_cell, [L, L, L])
        else:
            t.add_callback(fix_cell, [float(_) for _ in L])

    # Define slice.
    # We interpret --first N --last N as a request of step N
    if args.last == args.first and args.last is not None:
        args.last += 1
    sl = fractional_slice(args.first, args.last, args.skip, len(t))

    # Unfold if requested
    if args.unfold:
        tu = trajectory.Unfolded(t) #, fix_cm=True)
    else:
        tu = t

    # Fix CM and fold back
    if args.fix_cm:
        tu.add_callback(trajectory.fix_cm)
        tu.add_callback(trajectory.fold)

    # Here we could you a trajectory slice t[sl] but this will load
    # everything in ram (getitem doesnt provide a generator). This
    # will be fixed with python 3.
    ts = trajectory.Sliced(tu, sl)

    # Change number of particles
    if args.N > 0:
        def decimate_system(system, N):
            from atooms.system.particle import decimate
            system.particle = decimate(system.particle, N)
            return system
        ts.register_callback(decimate_system, args.N)
    # Change density and temperature
    if args.rho is not None:
        ts.register_callback(trajectory.decorators.set_density, args.rho)
    if args.temperature is not None:
        ts.register_callback(trajectory.decorators.set_temperature, args.temperature)
    # Change species layout if requested
    if args.species_layout is not None:
        ts.register_callback(trajectory.decorators.change_species, args.species_layout)
    # Sort by species id
    if args.species_sort:
        ts.register_callback(trajectory.decorators.sort)

    # We enforce regular periodicity; steps is None is trajectory is not periodic
    try:
        steps = check_block_size(ts.steps, ts.block_size, prune=True)
    except IndexError:
        print('# Warning: something wrong with periodicity check.')
        print('# We will proceed, but you should check the converted trajectory.')
        steps = ts.steps

    #
    # ---------------------
    # Trajectory conversion
    # ---------------------
    #
    include_list, exclude_list = [], []
    if len(args.fields_include) > 0:
        include_list = args.fields_include.split(',')
    if len(args.fields_exclude) > 0:
        exclude_list = args.fields_exclude.split(',')

    from atooms.trajectory import Trajectory
    if isinstance(out_class, str):
        out_class = Trajectory.formats[out_class]

    conv = ts.copy(fout=args.file_out, cls=out_class,
                   only=args.fields, exclude=exclude_list,
                   include=include_list, steps=steps)
    fout = conv.close()

    if args.ff:
        from atooms.trajectory.hdf5 import add_interaction_hdf5
        add_interaction_hdf5(fout, args.ff)

    if args.verbose and args.file_out != '/dev/stdout':
        print('# converted %s to %s' % (args.file_inp, fout))

    t.close()


def main_paste(args):
    """
    Correlate particles properties from trajectory files.

    Example:
    --------

    trj.py paste.py file1.xyz:radius file2.xyz.voronoi.xyz:volume
    """
    from atooms import trajectory as trj
    from atooms.core.utils import tipify

    f1, attr1 = args.file_inp[0].split(':')
    f2, attr2 = args.file_inp[1].split(':')
    if args.inp is None:
        fmt1, fmt2 = None, None
    else:
        fmt1, fmt2 = args.inp.split(',')
    t1 = trj.Trajectory(f1, fmt=fmt1)
    t2 = trj.Trajectory(f2, fmt=fmt2)

    # Define slice.
    # We interpret --first N --last N as a request of step N
    if args.last == args.first and args.last is not None:
        args.last += 1
    sl1 = fractional_slice(args.first, args.last, args.skip, len(t1))
    sl2 = fractional_slice(args.first, args.last, args.skip, len(t2))

    # Here we could you a trajectory slice t[sl] but this will load
    # everything in ram (getitem doesnt provide a generator). This
    # will be fixed with python 3.
    ts1 = trajectory.Sliced(t1, sl1)
    ts2 = trajectory.Sliced(t2, sl2)

    def array_fmt(arr):
        """Remove commas and [] from numpy array repr."""
        # Passing a scalar will trigger an error (gotcha: even
        # when casting numpy array to list, the elements remain of
        # numpy type and this function gets called! (4% slowdown)
        _fmt = '%g'
        try:
            return ' '.join([_fmt % x for x in arr])
        except:
            return _fmt % arr
            # except:
            #     return numpy.array2string(arr, precision=self.precision, separator=',')[1:-1]
    import numpy
    numpy.set_string_function(array_fmt, repr=False)

    for step, s1, s2 in trj.utils.paste(ts1, ts2):
        try:
            for i in range(len(s1.particle)):
                print(getattr(s1.particle[i], attr1),
                      getattr(s2.particle[i], attr2))
        except:
            print(getattr(s1, attr1), getattr(s2, attr2))

def scatter(args):
    """
    Write frames in trajectory to individual files of the same file format
    """
    from atooms import trajectory as trj
    from atooms.core.utils import mkdir

    for f in args.file_inp:
        fmt = args.inp
        fmt_out = fmt
        t = trj.Trajectory(f, fmt=fmt)

        # Define slice
        # We interpret --first N --last N as a request of step N
        if args.last == args.first and args.last is not None:
            args.last += 1
        sl = fractional_slice(args.first, args.last, args.skip, len(t))
        ts = trajectory.Sliced(t, sl)
        for i, system in enumerate(ts):
            f_out = args.file_out.format(step=ts.steps[i], frame=i,
                                         base=os.path.splitext(f)[0],
                                         ext=os.path.splitext(f)[1])
            mkdir(os.path.dirname(f_out))
            with trj.Trajectory(f_out, fmt=fmt_out, mode='w') as th_out:
                if args.fields is not None:
                    th_out.variables = args.fields.split(',')
                else:
                    th_out.variables = ts.variables
                th_out.write(system, step=ts.steps[i])


if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser(epilog=formats(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser = add_first_last_skip(parser)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')
    parser.set_defaults(func=lambda args: parser.print_help())
    subparsers = parser.add_subparsers()

    parser_convert = subparsers.add_parser('convert')
    parser_convert.add_argument(      '--fields', dest='fields', help='attributes fields')
    parser_convert.add_argument('-I', '--fields-include', dest='fields_include', type=str, default='', help='include patterns in fields')
    parser_convert.add_argument('-E', '--fields-exclude', dest='fields_exclude', type=str, default='', help='exclude patterns from fields')
    parser_convert.add_argument('-i', '--fmt-inp', dest='inp', help='input format')
    parser_convert.add_argument('-o', '--fmt-out', dest='out', help='output format for conversion')
    parser_convert.add_argument(      '--folder', dest='folder', action='store_true', help='force folder-based layout')
    parser_convert.add_argument('-F', '--ff', dest='ff', type=str, default='', help='force field file')
    parser_convert.add_argument(      '--flatten-steps',dest='flatten_steps', action='store_true', help='use sample index instead of steps')
    parser_convert.add_argument(      '--unfold',dest='unfold', action='store_true', help='unfold')
    parser_convert.add_argument(      '--fix-cm',dest='fix_cm', action='store_true', help='fix cm')
    parser_convert.add_argument(      '--side', dest='side', default=None, help='set cell side')
    parser_convert.add_argument(      '--density', dest='rho', type=float, default=None, help='new density')
    parser_convert.add_argument('-T', '--temperature', dest='temperature', type=float, default=None, help='new temperature')
    parser_convert.add_argument('-N', '--number-of-particles', dest='N', type=int, default=-1, help='new number of particles')
    parser_convert.add_argument(      '--precision', dest='precision', type=int, default=None, help='write precision')
    parser_convert.add_argument(      '--species',dest='species_layout', default=None, help='modify species layout (A, C, F)')
    parser_convert.add_argument(      '--sort-species',dest='species_sort', action='store_true', help='sort by species')
    parser_convert.add_argument(      '--seed', dest='seed', type=int, help='set seed of random number generator')
    parser_convert.add_argument(nargs=1, dest='file_inp', default='-', help='input file')
    parser_convert.add_argument(nargs='?', dest='file_out', default='-', help='output file')
    parser_convert.set_defaults(func=main)

    parser_info = subparsers.add_parser('info')
    parser_info.add_argument(      '--folder', dest='folder', action='store_true', help='force folder-based layout')
    parser_info.add_argument(      '--what', dest='fields', help='what info to show')
    parser_info.add_argument('-i', '--fmt-inp', dest='inp', help='input format')
    parser_info.add_argument(nargs='*', dest='file_inp', default='-', help='input file')
    parser_info.set_defaults(func=main_info)

    parser_paste = subparsers.add_parser('paste')
    parser_paste.add_argument('-i', '--fmt-inp', dest='inp', help='input format')
    parser_paste.add_argument(nargs=2, dest='file_inp', help='input files')
    parser_paste.set_defaults(func=main_paste)

    parser_scatter = subparsers.add_parser('scatter')
    parser_scatter.add_argument(      '--fields', dest='fields', help='attributes fields')
    parser_scatter.add_argument('-i', '--fmt-inp', dest='inp', help='input format')
    parser_scatter.add_argument('-o', '--output-file', dest='file_out', default='{base}_{frame}{ext}', help='output path (interpolated)')
    parser_scatter.add_argument(nargs='*', dest='file_inp', help='input file')
    parser_scatter.set_defaults(func=scatter)

    # parse argument lists
    args = parser.parse_args()
    args.func(args)
