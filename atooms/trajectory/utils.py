"""Useful functions to manipulate trajectories."""

import os
import tarfile
import numpy
import copy
from atooms.core.progress import progress


def gopen(filename, mode):
    """Open a file recognizing gzipped and bzipped files by extension."""
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode + 't')
    elif ext == '.bz2':
        import bz2
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)


def modify_fields(trajectory, fields=None, include=None, exclude=None):
    """
    Modify fields of a trajectory.

    Either provide a new list of fields, such as ['id', 'x', 'y'], or
    specify explicit patterns to exclude or include.
    """
    if fields is not None:
        # Reset the output format
        trajectory.fields = fields
    else:
        # Exclude and/or include lists of patterns from output format
        if exclude is not None:
            for pattern in exclude:
                if pattern in trajectory.fields:
                    trajectory.fields.remove(pattern)
        if include is not None:
            for pattern in include:
                if pattern not in trajectory.fields:
                    trajectory.fields.append(pattern)

    return trajectory


def convert(inp, out, fout, force=True, fields=None,
            exclude=None, include=None, steps=None):
    """
    Convert trajectory into a different format.

    `inp`: input trajectory object
    `out`: output trajectory class
    `fout`: output file

    If `out` is a string, we look for a matching trajectory format
    else we assume out is a trajectory class.
    If `out` is None, we rely on the factory guessing the format
    from the filename suffix.

    Return: name of converted trajectory file
    """
    from atooms.trajectory import Trajectory
    from atooms.trajectory.base import canonicalize_fields
    if isinstance(out, str):
        out_class = Trajectory.formats[out]
    else:
        out_class = out

    if fields is None and len(inp.fields) > 0 and include is None:
        # We automatically include all fields from the input trajectory
        # Since the output trajectory may have extra fields, we do should not overwrite them
        include = canonicalize_fields(inp.fields)

    if fout != '/dev/stdout' and (os.path.exists(fout) and not force):
        print('File exists, conversion skipped')
    else:
        # Make sure parent folder exists
        from atooms.core.utils import mkdir
        mkdir(os.path.dirname(fout))
        with out_class(fout, 'w') as conv:
            modify_fields(conv, fields, include, exclude)
            conv.precision = inp.precision
            conv.timestep = inp.timestep
            conv.block_size = inp.block_size
            # In python 3, zip returns a generator so this is ok
            #
            # for system, step in zip(inp, inp.steps):
            #     conv.write(system, step)
            #
            # In python 2, zipping t and t.steps will load everything
            # in RAM. In this case, it is better to use enumerate()
            if steps is None:
                for i, system in progress(enumerate(inp), total=len(inp)):
                    conv.write(system, inp.steps[i])
            else:
                # Only include requested steps (useful to prune
                # non-periodic trajectories)
                for step in steps:
                    idx = inp.steps.index(step)
                    conv.write(inp[idx], step)

    return fout


def split(inp, index='step', archive=False):
    """
    Split the trajectory into independent trajectory files, one per
    sample.
    """
    if archive:
        tar = tarfile.open(inp.filename + '.tar.gz', "w:gz")
    base, ext = os.path.splitext(inp.filename)

    for frame, step in enumerate(inp.steps):
        system = inp[frame]
        if index == 'step':
            filename = '%s-%09i%s' % (base, step, ext)
        elif index == 'frame' or index == 'sample':
            filename = '%s-%09i%s' % (base, frame, ext)
        else:
            raise ValueError('unknown option %s' % index)
        with inp.__class__(filename, 'w') as t:
            t.write(system, step)
        if archive:
            tar.add(filename)
            os.remove(filename)

    if archive:
        tar.close()


def get_block_size(data):
    """
    Return the size of the periodic block after which entries in
    `data` repeat. It is used to determine the block size in
    trajectories with logarithmic time spacing.
    """
    if len(data) < 2:
        return 1
    delta_old = 0
    delta_one = data[1] - data[0]
    iold = data[0]
    period = 1
    for ii in range(1, len(data)):
        i = data[ii]
        delta = i-iold
        # If we find that we repeat the increment between entries is
        # smaller than the previous iteration and it gets back to the
        # initial one (delat_one) then we found a block. We must
        # correct the +1 overshoot thus we subtract -1 to period
        if delta < delta_old and delta == delta_one:
            return period - 1
        else:
            period += 1
            iold = i
            delta_old = delta

    # We got to the end of the trajectory
    if len(data) != period:
        raise ValueError('something went wrong in block analysis')
    if data[1]-data[0] == data[-1]-data[-2]:
        # If the difference between steps is constant (euristically)
        # the period is one
        return 1
    else:
        # There is no periodicity, the block size is the whole trajectory
        return period


def check_block_size(steps, block_size, prune=False):
    """
    Perform some consistency checks on periodicity of non linear sampling.

    `block_size` is the number of frames composing a periodic block.
    If `prune` is True, the steps that do not match the first periodic
    block will be removed.

    Return a new list of steps that match the periodicity.

    Example:
    -------
    steps = [0, 1, 2, 4, 8, 9, 10, 12, 16]
    block_size = 4

    Note that in this case, len(steps) % block_size == 1, which is tolerated.
    """
    if block_size == 1:
        return None

    steps_local = copy.copy(steps)

    # Identify steps that do not match the first periodic block
    block = steps_local[0: block_size]
    ibl, jbl = 0, 0
    prune_me = []
    for i, step in enumerate(steps_local):
        step_expected = ibl * steps_local[block_size] + block[jbl]
        if step == step_expected:
            if jbl == block_size-1:
                # We are done with this block, we start over
                ibl += 1
                jbl = 0
            else:
                # We increment the index within the block
                jbl += 1
        else:
            prune_me.append(step)

    # Remove samples that do not conform with first block
    if prune and len(prune_me) > 0:
        print('#', len(prune_me), 'samples should be pruned')
        for step in prune_me:
            _ = steps_local.pop(steps_local.index(step))

    # Check if the number of steps is an integer multiple of
    # block period (we tolerate a rest of 1)
    rest = len(steps_local) % block_size
    if rest > 1:
        steps_local = steps_local[:-rest]
        print('# block was truncated')

    # Final test, after pruning spurious samples we should have a period
    # sampling, otherwise there was some error
    nbl = len(steps_local) // block_size
    for i in range(nbl):
        i0 = steps_local[i * block_size]
        current = steps_local[i * block_size: (i + 1) * block_size]
        current = [ii - i0 for ii in current]
        if not current == block:
            print('# periodicity issue at block %i out of %i' % (i, nbl))
            print('# current     :', current)
            print('# finger print:', block)
            raise IndexError('block does not match finger print')

    return steps_local


def dump(trajectory, what='pos'):
    """
    Dump coordinates as a list of (npart, ndim) numpy arrays if the
    trajectory is grandcanonical or as (nsteps, npart, ndim) numpy
    array if it is not grandcanonical.
    """
    if trajectory.grandcanonical:
        data = []
        for i, s in enumerate(trajectory):
            data[i].append(s.dump(what))
    else:
        data = numpy.zeros([len(trajectory.steps),
                            len(trajectory[0].particle),
                            len(trajectory[0].cell.side)])
        for i, s in enumerate(trajectory):
            data[i] = s.dump(what)

    return data


def field(trajectory, trajectory_field, field_name, frame, x_field=None):
    """
    Return the field specified by particle attribute `field_name` at a
    given `frame`.
    """
    if x_field is not None:
        raise DeprecationWarning('use field_name instead of x_field')
        field_name = x_field
    step = trajectory.steps[frame]
    try:
        index_field = trajectory_field.steps.index(step)
    except ValueError:
        return None
    x = []
    for pi in trajectory_field[index_field].particle:
        fi = getattr(pi, field_name)
        x.append(fi)
    return x


def paste(t1, t2):
    """
    Iterate simultaneously on two trajectories. Skip samples that
    exist in one trajectory and not in the other.

    Example:
    -------
    t1 = Trajectory(f1)
    t2 = Trajectory(f2)
    for s1, s2 in paste(t1, t2):
        pass
    """
    steps_1 = set(t1.steps)
    steps_2 = set(t2.steps)
    steps = sorted(steps_1 & steps_2)
    for step in steps:
        s1 = t1[t1.steps.index(step)]
        s2 = t2[t2.steps.index(step)]
        yield step, s1, s2


def is_cell_variable(trajectory, tests=1):
    """
    Simple test to check if cell changes.

    We compare the first frame to an integer number `tests` of other
    frames starting from the end of `trajectory`.
    """
    is_variable = False
    frames = len(trajectory)
    if tests > 0:
        skip = max(1, int(frames / float(tests)))
    else:
        skip = 1
    L0 = trajectory[0].cell.side
    for sample in range(frames-1, 0, -skip):
        L1 = trajectory[sample].cell.side
        if L0[0] != L1[0] or L0[1] != L1[1] or L0[2] != L1[2]:
            is_variable = True
            break
    return is_variable


def formats():
    """Return a string with the available trajectory formats."""
    from atooms import trajectory
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


def info(trajectory, keys=None):
    """Return a string with information about a `trajectory` instance."""
    from atooms.system.particle import distinct_species, composition

    if keys is None:
        # Default: full info
        txt = ''
        txt += 'path                 = %s\n' % trajectory.filename
        txt += 'format               = %s\n' % trajectory.__class__
        txt += 'frames               = %s\n' % len(trajectory)
        txt += 'megabytes            = %s\n' % (os.path.getsize(trajectory.filename) / 1e6)
        txt += 'particles            = %s\n' % len(trajectory[0].particle)
        txt += 'species              = %s\n' % ', '.join(distinct_species(trajectory[0].particle))
        txt += 'composition          = %s\n' % dict(composition(trajectory[0].particle))
        txt += 'density              = %s\n' % round(trajectory[0].density, 10)
        txt += 'cell side            = %s\n' % str(list(trajectory[0].cell.side))[1: -1]
        txt += 'cell volume          = %s\n' % trajectory[0].cell.volume
        if len(trajectory) > 1:
            txt += 'steps                = %s\n' % trajectory.steps[-1]
            txt += 'duration             = %s\n' % trajectory.times[-1]
            txt += 'timestep             = %s\n' % trajectory.timestep
            txt += 'block size           = %s\n' % trajectory.block_size
            if trajectory.block_size == 1:
                txt += 'steps between frames = %s\n' % (trajectory.steps[1]-trajectory.steps[0])
                txt += 'time between frames  = %s\n' % (trajectory.times[1]-trajectory.times[0])
            else:
                txt += 'block steps          = %s\n' % trajectory.steps[trajectory.block_size-1]
                txt += 'block                = %s\n' % ([trajectory.steps[i] for i in range(trajectory.block_size)])
            txt += 'grandcanonical       = %s' % trajectory.grandcanonical
        return txt

    else:
        # Selected infos.
        # TODO: of course, it would be cleaner to have a little class for that
        outs = []
        for key in keys.split(','):
            if key == 'path':
                outs.append(trajectory.filename)
            elif key == 'format':
                outs.append(trajectory.__class__)
            elif key == 'frames':
                outs.append(len(trajectory))
            elif key == 'megabytes':
                outs.append(os.path.getsize(trajectory.filename) / 1e6)
            elif key == 'particles':
                outs.append(len(trajectory[0].particle))
            elif key == 'species':
                outs.append(', '.join(distinct_species(trajectory[0].particle)))
            elif key == 'composition':
                outs.append(dict(composition(trajectory[0].particle)))
            elif key == 'cell density':
                outs.append(round(trajectory[0].density, 10))
            elif key == 'cell side':
                outs.append(str(list(trajectory[0].cell.side))[1: -1])
            elif key == 'cell volume':
                outs.append(trajectory[0].cell.volume)
            elif key == 'steps':
                outs.append(trajectory.steps[-1])
            elif key == 'duration':
                outs.append(trajectory.times[-1])
            elif key == 'timestep':
                outs.append(trajectory.timestep)
            elif key == 'block size':
                outs.append(trajectory.block_size)
            elif key == 'steps between frames':
                outs.append(trajectory.steps[1]-trajectory.steps[0])
            elif key == 'time between frames':
                outs.append(trajectory.times[1]-trajectory.times[0])
            elif key == 'block steps':
                outs.append(trajectory.steps[trajectory.block_size-1])
            elif key == 'block':
                outs.append([trajectory.steps[i] for i in range(trajectory.block_size)])
            elif key == 'grandcanonical':
                outs.append(trajectory.grandcanonical)

        txt = ''
        fmt = '%%-%ds : %%s\n' % (max([len(key) for key in keys.split(',')]))
        for key, out in zip(keys.split(','), outs):
            txt += fmt % (key, out)
        
        return txt.strip('\n')
