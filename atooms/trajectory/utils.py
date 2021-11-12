"""Useful functions to manipulate trajectories."""

import os
import copy
import tarfile
import numpy
import warnings

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


def file_index(fh, size=None):
    """Lightweight file indexing for trajectories via tell/seek"""
    header = []
    block = []
    block_size = []
    if size is None:
        def size(fh, data, line):
            npart = int(data)
            return 1, npart

    # Make sure we start from the beginning
    fh.seek(0)
    while True:
        line = fh.tell()
        data = fh.readline()
        # We break if file is over or we found an empty line
        if not data:
            break

        # Get size of header and block
        try:
            header_size, this_block_size = size(fh, data, line)
            header.append(line)
        except ValueError:
            raise IOError('malformed file [{}]'.format(fh.filename))

        # Skip header_size lines (if zero none will be skipped)
        for _ in range(header_size):
            fh.readline()

        # Skip block_size lines, making sure we have
        # read precisely that number of lines
        line = fh.tell()
        for _ in range(this_block_size):
            data = fh.readline()

        # Store first line /after/ we have read the frame
        # making sure the last we read was not emtpy
        # Note that readline() returns an empty string on EOF
        if len(data) > 0:
            block.append(line)
            block_size.append(this_block_size)
        else:
            raise IOError('malformed file [%s]', fh.filename)

    # TODO: leave fh as it was
    return header, block, block_size


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
    # Linear sampling
    if block_size == 1:
        return None

    # This is a single non-linear block
    if block_size == len(steps):
        return None

    steps_local = copy.copy(steps)

    # Identify steps that do not match the first periodic block
    block = steps_local[0: block_size]
    ibl, jbl = 0, 0
    prune_me = []
    for _, step in enumerate(steps_local):
        offset = block[0] if ibl > 0 else 0
        step_expected = ibl * (steps_local[block_size] - offset) + block[jbl]
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
        for step in prune_me:
            _ = steps_local.pop(steps_local.index(step))

    # Check if the number of steps is an integer multiple of
    # block period (we tolerate a rest of 1)
    rest = len(steps_local) % block_size
    if rest > 1:
        steps_local = steps_local[:-rest]
        warnings.warn('truncated block')

    # Final test, after pruning spurious samples we should have a period
    # sampling, otherwise there was some error
    nbl = len(steps_local) // block_size
    for i in range(nbl):
        # We test that the difference between the finger print and the
        # first sample in the block is constant
        diff_last = None
        for j in range(len(block)):
            diff = steps_local[i*block_size + j] - block[j]
            if diff_last is None:
                diff_last = diff
            if diff_last != diff:
                raise IndexError('block does not match finger print {}'.format(block))
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
    for sample in range(frames-1, -1, -skip):
        L1 = trajectory[sample].cell.side
        if numpy.any(L0 != L1):
            is_variable = True
            break
    return is_variable


def is_semigrandcanonical(trajectory, tests=1):
    """
    Simple test to check if a trajectory is semigrandcanonical.
    i.e. if the chemical concentrations fluctuate.

    We compare the first frame to an integer number `tests` of other
    frames starting from the end of `trajectory`.
    """
    # This is adapted from is_cell_variable()
    is_variable = False
    if tests > 0:
        skip = max(1, int(len(trajectory) / float(tests)))
    else:
        skip = 1
    from atooms.system.particle import composition
    x0 = composition(trajectory[0].particle)
    for sample in range(len(trajectory)-1, -1, -skip):
        x1 = composition(trajectory[sample].particle)
        is_variable = False
        for sp in x0:
            if x0[sp] != x1[sp]:
                is_variable = True
                break
        if is_variable:
            break
    return is_variable


def is_grandcanonical(trajectory, tests=1):
    """
    Simple test to check if a trajectory is grandcanonical.
    i.e. if the number of particles fluctuates.

    We compare the first frame to an integer number `tests` of other
    frames starting from the end of `trajectory`.
    """
    # This is adapted from is_cell_variable()
    # and basically the same code as is_semigrandcanonical()
    is_variable = False
    if tests > 0:
        skip = max(1, int(len(trajectory) / float(tests)))
    else:
        skip = 1
    N0 = len(trajectory[0].particle)
    for sample in range(len(trajectory)-1, 0, -skip):
        N1 = len(trajectory[sample].particle)
        is_variable = False
        if N0 != N1:
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
    system = trajectory[0]
    if keys is None:

        # Default: full info
        txt = ''
        txt += 'path                 = %s\n' % trajectory.filename
        txt += 'format               = %s\n' % trajectory.__class__
        txt += 'frames               = %s\n' % len(trajectory)
        txt += 'megabytes            = %s\n' % (os.path.getsize(trajectory.filename) / 1e6)
        txt += 'particles            = %s\n' % len(system.particle)
        txt += 'species              = %s\n' % ', '.join(distinct_species(system.particle))
        txt += 'composition          = %s\n' % dict(composition(system.particle))
        txt += 'size dispersion      = {}\n'.format((numpy.std([p.radius for p in system.particle]) / numpy.mean([p.radius for p in system.particle])))
        txt += 'density              = %s\n' % round(system.density, 10)
        if system.cell is not None:
            txt += 'cell side            = %s\n' % str(list(system.cell.side))[1: -1]
            txt += 'cell volume          = %s\n' % system.cell.volume
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
                outs.append(len(system.particle))
            elif key == 'species':
                outs.append(', '.join(distinct_species(system.particle)))
            elif key == 'composition':
                outs.append(dict(composition(system.particle)))
            elif key == 'cell density':
                outs.append(round(system.density, 10))
            elif key == 'cell side':
                outs.append(str(list(system.cell.side))[1: -1])
            elif key == 'cell volume':
                outs.append(system.cell.volume)
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
