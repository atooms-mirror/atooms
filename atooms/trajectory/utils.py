
"""Useful functions to manipulate trajectories."""

import os
import tarfile
import gzip
import re

def gopen(filename, mode):
    """Open a file recognizing gzipped files by extension."""
    import gzip
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        return gzip.open(filename, mode)
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

def convert(inp, out, fout='', tag='', prefix='', force=True, fmt=None, exclude=[], include=[], stdout=False, callback=None, args={}):
    # TODO: check these dangerous defaults
    """Convert trajectory into a different format.

    inp: input trajectory object
    out: output trajectory class
    tag: optional string to be prepended before the output suffix

    Return: name of converted trajectory file
    """
    # TODO: convert metadata (interaction etc) !
    # If the input trajectory lies in a directory, the new trajectory is located
    # in a companion directory prefixed by tag. The basename is config

    # # Check that we have some files there
    # if len(inp) == 0:
    #     raise IOError('no samples in trajectory (%s)' % inp.filename)

    if stdout:
        filename = '/dev/stdout'
    else:
        if len(fout) > 0:
            filename = fout
        elif os.path.isdir(inp.filename):
            # if tag == '':
            #     if prefix == '':
            #         tag = '-conv'
            d = os.path.dirname(inp.filename)
            b = os.path.basename(inp.filename)
            dirname = os.path.join(d, prefix + b) + tag
            filename = dirname + '.' + out.suffix
            # from pyutils.utils import mkdir
            # mkdir(dirname)
            # filename = dirname + '/config.' + out.suffix
        else:
            filename = os.path.splitext(inp.filename)[0] + tag + '.' + out.suffix    

    if not os.path.exists(filename) or force:
        with out(filename, 'w') as conv:
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

    return filename


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
    period = 0
    for ii in range(1, len(data)):
        i = data[ii]
        delta = i-iold
        if delta < delta_old and delta == delta_one and abs(delta-delta_old)>2:
            return period
        else:
            period += 1
            iold = i
            delta_old = delta
    return 1
