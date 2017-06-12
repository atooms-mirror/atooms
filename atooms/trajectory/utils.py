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
                if not pattern in trj.fmt:
                    trj.fmt.append(pattern)
                    
    return trj

def convert(inp, out, fout, tag='', prefix='', force=True, fmt=None,
            exclude=None, include=None, callback=None, args=None):
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
            print fmt, conv.fmt, conv._fmt, include, exclude
            format_output(conv, fmt, include, exclude)
            conv.precision = inp.precision
            conv.timestep = inp.timestep
            conv.block_size = inp.block_size
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
    file_steps.sort(key=lambda a: a[1])
    new_files = [a[0] for a in file_steps]
    new_steps = [a[1] for a in file_steps]
    return new_files, new_steps

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

def check_block_size(steps, block_size):
    """
    Perform some consistency checks on periodicity of non linear sampling.

    `block_size` is the number of frames composing a periodic block.

    Example:
    -------
    steps = [0, 1, 2, 4, 8, 9, 10, 12, 16]
    block_size = 4

    Note that in this case, len(steps) % block_size == 1, which is tolerated.
    """
    if block_size == 1:
        return

    block = steps[0:block_size]
    ibl = 0
    jbl = 0
    prune_me = []
    for i in steps:
        j = ibl * steps[block_size] + block[jbl]
        if i == j:
            if jbl == block_size-1:
                # We are done with this block, we start over
                ibl += 1
                jbl = 0
            else:
                # We increment the index within the block
                jbl += 1
        else:
            prune_me.append(i)

    if len(prune_me) > 0:
        print '\n# ', len(prune_me), ' samples will be pruned'

    for p in prune_me:
        pp = steps.index(p)

    # check if the number of steps is an integer multiple of
    # block period (we tolerate a rest of 1)
    rest = len(steps) % block_size
    if rest > 1:
        steps = steps[:-rest]
        print 'block was truncated'

    # final test, after pruning spurious samples we should have a period
    # sampling, otherwise there was some error
    nbl = len(steps) / block_size
    for i in range(nbl):
        i0 = steps[i*block_size]
        current = steps[i*block_size:(i+1)*block_size]
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
        data = numpy.zeros([len(trajectory.steps),
                            len(trajectory[0].particle),
                            len(trajectory[0].cell.side)])
        for i, s in enumerate(trajectory):
            data[i] = s.dump(what)

    return data

def field(trajectory, trajectory_field, field, sample):
    step = trajectory.steps[sample]
    try:
        index_field = trajectory_field.steps.index(step)
    except:
        return None
    x = []
    for pi in trajectory_field[index_field].particle:
        fi = getattr(pi, field)
        #fi = [float(x) for x in fi.split(',')]
        x.append(fi)
    return x

def tzip(t1, t2):
    """
    Iterate simultaneously on two trajectories. Skip samples that
    exist in one trajectory and not in the other.

    Example:
    -------
    t1 = Trajectory(f1)
    t2 = Trajectory(f2)
    for s1, s2 in tzip(t1, t2):
        pass
    """
    steps_1 = set(t1.steps)
    steps_2 = set(t2.steps)
    steps = steps_1 & steps_2
    for step in steps:
        s1 = t1[t1.steps.index(step)]
        s2 = t2[t2.steps.index(step)]
        yield step, s1, s2

def dump(*fileinp):

    import sys
    from atooms import trajectory as trj                          

    particle=False
    t_list, attr_list = [], []
    for f in fileinp:
        f_i, attr_i = f.split(':')
        if 'particle' in attr_i:
            attri_i = attr_i.split('particle.')[1]
            particle = True
        t_list.append(Trajectory(f_i))
        attr_list.append(attr_i)

    steps = None
    for t in t_list:
        steps_i = set(t.steps)
        if steps is None:
            steps = steps_i
        else:
            steps = steps & steps_i
    
    for step in steps:
        if not particle:
            txt = ''
            for t, attr in zip(t_list, attr_list):
                s = t[t.steps.index(step)]
                txt += str(getattr(s, attr))
            print txt

def available_formats():
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
            txt += 'block steps          = %s\n' % trajectory.steps[trajectory.block_size-1]
            txt += 'block                = %s\n' % ([trajectory.steps[i] for i in range(trajectory.block_size)])
        txt += 'grandcanonical       = %s' % trajectory.grandcanonical
    print txt

def read(file_inp, inp=None):
    from atooms.utils import Timer

    t = trajectory.Trajectory(file_inp, fmt=inp)
    t = Timer()
    t.start()
    for s in trajectory:
        pass
    t.stop()
    print t.wall_time, t.wall_time / 60.

def main(file_inp, file_out, inp=None, out=None, folder=False,
         precision=None, seed=None, side=None, rho=None,
         temperature=None, alphabetic_ids=None, tag='', fmt=None,
         fmt_include='', fmt_exclude='', ff=None, first=None,
         last=None, skip=1):
    """Convert trajectory `file_inp` to `file_out`."""
    import random
    from atooms import trajectory
    from atooms.utils import fractional_slice, add_first_last_skip

    if file_out == '-':
        file_out = '/dev/stdout'

    if folder:
        t = trajectory.folder.Foldered(file_inp, cls=inp)
    else:
        t = trajectory.Trajectory(file_inp, fmt=inp)

    # If no output format is provided we use the input one
    if out is None:
        out_class = t.__class__
    else:
        out_class = out

    if precision is not None:
        t.precision = precision

    if flatten_steps:
        t.steps = range(1,len(t)+1)

    # Reset random number generator
    if seed:
        random.seed(seed)

    # Trick to allow some trajectory formats to set the box side.
    # This way the cell is defined as we read the sample (callbacks
    # will not do that).
    if side is not None:
        t._side = side

    # Define slice.
    # We interpret --first N --last N as a request of step N
    if last == first and last is not None:
        last += 1
    sl = fractional_slice(first, last, skip, len(t))
    # Here we could you a trajectory slice t[sl] but this will load
    # everything in ram (getitem doesnt provide a generator). This
    # will be fixed with python 3.
    ts = trajectory.Sliced(t, sl)

    # Change density and temperature
    if rho is not None:
        ts.register_callback(trajectory.decorators.set_density, rho)
    if temperature is not None:
        ts.register_callback(trajectory.decorators.set_temperature, temperature)

    # We always normalize species id's using fortran convention
    ts.register_callback(trajectory.decorators.normalize_id, alphabetic_ids)

    # Trajectory conversion
    fout = trajectory.convert(ts, out_class, file_out,
                              tag=tag, fmt=fmt,
                              include=fmt_include.split(','),
                              exclude=fmt_exclude.split(','))

    if ff:
        from atooms.trajectory.hdf5 import add_interaction_hdf5
        add_interaction_hdf5(fout, ff)
    
    if file_out != '/dev/stdout':
        print '%s' % fout

    t.close()
