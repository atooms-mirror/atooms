"""Useful functions to manipulate trajectories."""

import os
import tarfile
import gzip
import re
import numpy

def gopen(filename, mode):
    """Open a file recognizing gzipped and bzipped files by extension."""
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode)
    elif ext == '.bz2':
        import bz2
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)

def format_output(trj, fmt=None, include=None, exclude=None):
    """
    Modify output format of an input trajectory.

    Either provide a new format, such as ['id', 'x', 'y'], or
    specify explicit patterns to exclude or include.
    """
    if fmt is not None:
        # Reset the output format
        trj.fmt = fmt
    else:
        # Exclude and/or include lists of patterns from output format
        if exclude is not None:
            for pattern in exclude:
                if pattern in trj.fmt:
                    trj.fmt.remove(pattern)
        if include is not None:
            for pattern in include:
                if pattern in trj.fmt:
                    trj.fmt.append(pattern)

    return trj

def convert(inp, out, fout, tag='', prefix='', force=True, fmt=None, exclude=[], include=[], callback=None, args={}):
    # TODO: check these dangerous defaults
    """Convert trajectory into a different format.

    inp: input trajectory object
    out: output trajectory class
    fout: output file

    Return: name of converted trajectory file
    """
    # TODO: convert metadata (interaction etc) !
    # If the input trajectory lies in a directory, the new trajectory is located
    # in a companion directory prefixed by tag. The basename is config

    # If out is a string, we look for a matching trajectory format
    # else we assume out is a trajectory class.
    # If out is None, we rely on the factory guessing the format
    # from the filename suffix.
    from atooms.trajectory import Trajectory
    if isinstance(out, basestring):
        out_class = Trajectory.formats[out]
    else:
        out_class = out

    if fout != '/dev/stdout' and (os.path.exists(fout) and not force):
        print 'File exists, conversion skipped'
    else:
        with out_class(fout, 'w') as conv:
            format_output(conv, fmt, include, exclude)
            conv.precision = inp.precision
            conv.timestep = inp.timestep
            conv.block_period = inp.block_period
            # TODO: Zipping t, t.steps is causing a massive mem leak!
            # In python <3 zip returns a list, not a generator! Therefore this
            # for system, step in zip(inp, inp.steps):
            #     conv.write(system, step)
            # will use a lot of RAM! Workarounds (in order of personal preference)
            # 1. zip is a generator in python 3
            # 2. use enumerate instead and grab the step from inp.steps[i]
            # 3. add an attribute system.step for convenience
            for i, system in enumerate(inp):
                if callback is not None:
                    system = callback(system, args)
                conv.write(system, inp.steps[i])

    return fout

def split(inp, selection=slice(None), index='step', archive=False):
    """Split the trajectory into independent trajectory files, one per sample."""
    if archive:
        tar = tarfile.open(inp.filename + '.tar.gz', "w:gz")
    base, ext = os.path.splitext(inp.filename)

    # TODO: fix zipping of steps
    for system, step, sample in zip(inp, inp.steps, inp.samples):
        if index == 'step':
            filename = '%s-%09i%s' % (base, step, ext)
        elif index == 'sample':
            filename = '%s-%09i%s' % (base, sample, ext)
        else:
            raise ValueError('unknown option %s' % index)
        with inp.__class__(filename, 'w') as t:
            t.write(system, step)
        if archive:
            tar.add(filename)
            os.remove(filename)

    if archive:
        tar.close()

def get_step_from_path(f):
    s = re.search(r'%s-(\d*)' % basename, f)
    if s:
        return int(s.group(1))

def sort_files_steps(files, steps):
    file_steps = zip(files, steps)
    file_steps.sort(key = lambda a : a[1])
    new_files = [a[0] for a in file_steps]
    new_steps = [a[1] for a in file_steps]
    return new_files, new_steps

def get_period(data):
    if len(data) < 2:
        return 1
    delta_old = 0
    delta_one = data[1] - data[0]
    iold = data[0]
    period = 1
    for ii in range(1, len(data)):
        i = data[ii]
        delta = i-iold
        if delta < delta_old and delta == delta_one: # and abs(delta-delta_old)>2:
            return period
        else:
            period += 1
            iold = i
            delta_old = delta
    if data[1]-data[0] == data[-1]-data[-2]:
        return 1
    else:
        return period

def check_block_period(steps, block_period):
    """Perform some consistency checks on periodicity of non linear sampling."""
    if block_period == 1:
        return

    block = steps[0:block_period]
    ibl = 0
    jbl = 0
    prune_me = []
    for i in steps:
        j = ibl*steps[block_period-1] + block[jbl]
        if i == j:
            jbl += 1
            if jbl == block_period-1:
                ibl += 1
                jbl = 0
        else:
            prune_me.append(i)

    if len(prune_me) > 0:
        print '\n# ', len(prune_me), ' samples will be pruned'

    for p in prune_me:
        pp = steps.index(p)

    # check if the number of steps is an integer multiple of
    # block period (we tolerate a rest of 1)
    rest = len(steps) % block_period
    if rest > 1:
        steps = steps[:-rest]
        #raise ValueError('block was truncated')

    # final test, after pruning spurious samples we should have a period
    # sampling, otherwise there was some error
    nbl = len(steps) / block_period
    for i in range(nbl):
        i0 = steps[i*block_period]
        current = steps[i*block_period:(i+1)*block_period]
        current = [ii-i0 for ii in current]
        if not current == block:
            print 'periodicity issue at block %i out of %i' % (i, nbl)
            print 'current     :', current
            print 'finger print:', block
            raise ValueError('block does not match finger print')


def dump(trajectory, what='pos'):
    """Dump coordinates as a list of (npart, ndim) numpy arrays if the
    trajectory is grandcanonical or as (nsteps, npart, ndim) numpy
    array if it is not grandcanonical.
    """
    if trajectory.grandcanonical:
        data = []
        for i, s in enumerate(trajectory):
            data[i].append(s.dump(what))
    else:
        data = numpy.zeros([len(trajectory.steps), \
                            len(trajectory[0].particle), \
                            len(trajectory[0].cell.side)])
        for i, s in enumerate(trajectory):
            data[i] = s.dump(what)

    return data
