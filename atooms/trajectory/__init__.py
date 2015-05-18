# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import numpy
import tarfile
import logging
import warnings

from decorators import *

class TrajectoryBase(object):

    """Trajectory base class"""

    suffix = None

    def __init__(self, filename, mode='r'):
        """When mode is 'r', it must set the list of available steps."""
        self.filename = filename
        self.mode  = mode
        self.steps = []

        # These are cached properties
        self._grandcanonical = None
        self._timestep = None
        self._block_period = None

        self._initialized_write = False
        self._initialized_read = False

        # These two are needed by unfold, system might be dropped
        self._system = None
        self._pos_unf = None

    def __len__(self):
        return len(self.steps)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __iter__(self):
        for i in range(len(self.steps)):
            yield self.read(i)

    def __getitem__(self, key):
        if isinstance(key, slice):
            raise TypeError("Not ready for slicing yet, use Sliced decorator instead")

        elif isinstance(key, int):
            if key < 0:
                key += len(self)
            if key >= len(self):
                raise IndexError("Index (%d) is out of range." % key)
            return self.read(key)

        else:
            raise TypeError("Invalid argument type [%s]" % type(key))

    def close(self):
        pass

    @property
    def samples(self):
        warnings.warn('iterate instead of using samples') #, DeprecationWarning)
        return range(len(self.steps))

    # Read and write implement the following template.
    #
    # 1. read_init() and write_init() are called only once to
    # intialize data structures or grab metadata (stuff that doesn't
    # change) 
    #
    # 2. read_sample() and write_sample() are used to actually
    # read/write a system
    #
    # Additionally, write_sample append the stp to the step list
    #
    # In future implementation, we might pass a list of objects to be
    # written, to store for instance, integrator data and so on.

    def read(self, index):
        if not self._initialized_read:
            self.read_init()
            self._initialized_read = True
        return self.read_sample(index)

    def write(self, system, step, ignore=[]):
        if self.mode == 'r':
            raise IOError('trajectory file not open for writing')
        if not self._initialized_write:
            self.write_init(system)
            self._initialized_write = True
        self.write_sample(system, step, ignore)
        # Step is added last, sample index starts from 0 by default
        self.steps.append(step)

    def read_init(self):
        """It may setup data structures need by the trajectory. Need not be implemented."""
        pass

    def write_init(self, system):
        """It may setup data structures need by the trajectory. Need not be implemented."""
        pass

    # These methods must be implemented by subclasses
    def read_sample(self, index): 
        """It must return the sample (system) with the given index"""
        raise NotImplementedError()
        
    def write_sample(self, system, step, ignore=[]):
        """It must write a sample (system) to disk. Noting to return."""
        raise NotImplementedError()

    # To read/write timestep and block period sublcasses may implement
    # these methods. The default is dt=1 and blockperiod determined dynamically.
    def read_timestep(self): 
        return 1.0

    def write_timestep(self, value): 
        pass

    def read_blockperiod(self): 
        return None

    def write_blockperiod(self, value): 
        pass

    @property
    def timestep(self):
        if self._timestep is None:
            self._timestep = self.read_timestep()
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self.write_timestep(value)
        self._timestep = value

    @property
    def block_period(self):
        if self._block_period is None:
            self._block_period = self.read_blockperiod()
        if self._block_period is None:
            # If period is still None (read_blockperiod is not
            # implemented) we determine it dynamically
            self._block_period = get_period(self.steps)
        return self._block_period

    @block_period.setter
    def block_period(self, value):
        self._block_period = value
        self.write_blockperiod(value)

    def _check_block_period(self):
        """Perform some consistency checks on periodicity of non linear sampling."""
        if self.block_period == 1:
            return
        block = self.steps[0:self.block_period]

        ibl = 0
        jbl = 0
        prune_me = []
        for k, i in enumerate(self.steps):
            j = ibl*self.steps[self.block_period] + block[jbl]
            if i == j:
                jbl += 1
                if jbl == self.block_period:
                    ibl += 1
                    jbl = 0
            else:
                prune_me.append(i)

        if len(prune_me) > 0:
            print '\n# ', len(prune_me), ' samples will be pruned'

        for p in prune_me:
            pp = self.steps.index(p)
            a = self.steps.pop(pp)

        # check if the number of steps is an integer multiple of
        # block period (we tolerate a rest of 1)
        rest = len(self.steps) % self.block_period
        if rest > 1:
            self.steps = self.steps[:-rest]
            #raise ValueError('block was truncated')

        # final test, after pruning spurious samples we should have a period
        # sampling, otherwise there was some error
        nbl = len(self.steps) / self.block_period
        for i in range(nbl):
            i0 = self.steps[i*self.block_period]
            current = self.steps[i*self.block_period:(i+1)*self.block_period]
            current = [ii-i0 for ii in current]
            if not current == block:
                print 'periodicity issue at block %i out of %i' % (i, nbl)
                print 'current     :', current
                print 'finger print:', block
                raise ValueError('block does not match finger print')

    # Some additional useful properties

    @property
    def grandcanonical(self): 
        # In subclasses, cache it for efficiency, since we might have to discover it
        if self._grandcanonical is None:
            self._grandcanonical = False
        return self._grandcanonical

    @property
    def times(self):
        """All available times."""
        return [s*self.timestep for s in self.steps]

    @property
    def time_total(self):
        """Total simulation time."""
        return self.steps[-1] * self.timestep

    def time_when_msd_is(self, msd_target, sigma=1.0):
        """Estimate the time when the MSD reaches target_msd in units of sigma^2.
        Bounded by the actual maximum time of trajectory tmax.
        """
        raise NotImplementedError('time_when_msd_is is broken')
        # TODO: when using decorators, this will introduce a diamond, let's move it to a generic function in post processing 
        # self._unfold()
        # msd_total = numpy.sum((self._pos_unf[-1] - self._pos_unf[0])**2) / self._pos_unf[0].shape[0]
        # return min(1.0, msd_target * sigma**2 / msd_total) * self.steps[-1] * self.timestep

    # def _unfold(self):
    #     # unfold should be handled only via the decorator.
    #     # We should allow subclasses to use internal data to avoid explicit unfolding.
    #     # Like cell image data or unfolded coordinates dumps
    #     if self._pos_unf:
    #         return
    #     self._system = self.read_initial_state()
    #     # do not include initial state, is it ok?
    #     pos = [] #[numpy.array([p.position for p in self._system.particle])]
    #     for sy in self:
    #         pos.append(numpy.array([p.position for p in sy.particle]))
    #         #pos[i] numpy.array([p.position for p in sy.particle]))
    #    # pos = [numpy.ndarray([p.position for p in self.read_sample(i).particle]) for i in self.samples]
    #     self._pos_unf = {}
    #     for i, p in zip(self.samples, _pbc_unfold(pos, self._system.cell.side)):
    #         self._pos_unf[i] = p


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

# Useful functions to manipulate trajectories

def convert(inp, out, tag='', ignore=[]):
    """Convert trajectory into a different format.

    inp: input trajectory object
    out: output trajectory class
    tag: optional string to be prepended before the output suffix

    Return: name of converted trajectory file
    """
    # TODO: strip trailing slash !
    # TODO: convert metadata (interaction etc) !
    filename = os.path.splitext(inp.filename)[0] + tag + '.' + out.suffix
    with out(filename, 'w') as conv:
        conv.timestep = inp.timestep
        conv.block_period = inp.block_period
        # In python <3 zip returns a list, not a generator! Therefore this
        # for system, step in zip(inp, inp.steps):
        #     conv.write(system, step, ignore=ignore)
        # will use a lot of RAM! Workarounds (in order of personal preference)
        # 1. zip is a generator in python 3
        # 2. use enumerate instead and grab the step from inp.steps[i]
        # 3. add an attribute system.step for convenience
        for i, system in enumerate(inp):
            conv.write(system, inp.steps[i], ignore=ignore)

    return filename

def split(inp, selection=slice(None), index='step', archive=False):
    """Split the trajectory into independent trajectory files, one per sample."""
    if archive:
        tar = tarfile.open(inp.filename + '.tar.gz', "w:gz")
    base, ext = os.path.splitext(inp.filename)

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

class SuperTrajectory(TrajectoryBase):

    """ Collection of subtrajectories """

    # The approach is inefficient for large number of
    # subtrajectories (init takes some time)
    # TODO; try not to keep all files open at a time !

    def __init__(self, subtrajectories, mode='r', timestep=1.0, variable=False, periodic=True):
        f = os.path.dirname(subtrajectories[0].filename) + '/'
        super(SuperTrajectory, self).__init__(f, mode)
        # By default we get info from the first trajectory in list
        #self.trajectory = subtrajectories[0]
        self.subtrajectories = subtrajectories
        self._timestep = timestep
        self.variable = variable

        # This is a bad hack to tolerate rumd hydiosincracies
        # Enforce all subt have the same length, the last one gets dropped
        if not self.variable:
            l = len(self.subtrajectories[0])
            prune = []
            for i in range(len(self.subtrajectories)):
                if l != len(self.subtrajectories[i]):
                    prune.append(i)
            for i in prune:
                self.subtrajectories.pop(i)

        # Make sure subtrajectories are sorted by increasing step
        self.subtrajectories.sort(key=lambda x : x.steps[0])

        # Mapping between samples and files that store them.
        # In simple cases we could get away with a modulo but this would work
        # e.g. when subtrajectories have different length
        self.steps = []
        self._map = []
        for i, t in enumerate(self.subtrajectories):
            # Check that neighboring subtrajectories
            # do not overlap at first/last step
            # The last block of course is not checked.
            if i > 0:
                if t.steps[0] == self.steps[-1]:
                    self.steps.pop(-1)
                    self._map.pop(-1)

            for sample, step in zip(t.samples, t.steps):
                self._map.append((i, sample))
                self.steps.append(step)

        # TODO: fix last block getting trimmed by check_block_period with rumd

        # Now it's better check the block period
        # Only check that if sampling is expected to be periodic.
        if periodic:
            self._check_block_period()

    def read_sample(self, sample):
        i, j = self._map[sample]
        return self.subtrajectories[i].read_sample(j)


# Factory method which mimics an abstract factory class
__factory_map = {}

from xyz import TrajectoryXYZ, TrajectoryPDB, TrajectoryNeighbors, TrajectoryXYZIkeda2
__factory_map['xyz'] = TrajectoryXYZ
__factory_map['pdb'] = TrajectoryPDB

from hoomd import TrajectoryHOOMD
__factory_map['tgz'] = TrajectoryHOOMD

from rumd import TrajectoryRUMD, SuperTrajectoryRUMD
from lammps import TrajectoryLAMMPS

try:
    from hdf5 import TrajectoryHDF5
    __factory_map['h5'] = TrajectoryHDF5
    __factory_map['dat'] = TrajectoryHDF5
except:
    pass

# TODO: trajectories should implement a method to check if a file
# is of their own format or not, to avoid relying on suffix
# check out http://stackoverflow.com/questions/456672/class-factory-in-python

# def Trajectory(filename, fmt=None, **kwargs):
# Make interface consistent with base, fmt is specified thorugh extension only
# It should be called Trajectory(fname, mode='w', fmt='h5') or
# Trajectory(fname, 'w', fmt='h5')
def Trajectory(filename, mode='r', fmt='' ): #, *args, **kwargs):
    """ Factory class shortcut """
    suffix = os.path.splitext(filename)[-1].replace('.', '')
    if not suffix in __factory_map:
        # always try hdf5 to accomodate non standard suffixes
        # we could even search within the path... :-)
        suffix = 'h5'
        #raise ValueError('unknown file type %s' % suffix)
        
    return __factory_map.get(suffix)(filename, mode)
#    return __factory_map.get(suffix)(filename, **kwargs)
