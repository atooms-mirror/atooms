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
            raise TypeError("Invalid argument type")

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
        return self.read_sample(index)

    def write(self, system, step):
        if self.mode == 'r':
            raise IOError('trajectory file not open for writing')
        if not self._initialized_write:
            self.write_init(system)
        self.write_sample(system, step)
        # Step is added last, sample index starts from 0 by default
        self.steps.append(step)

    # These methods must be implemented by subclasses
    def read_init(self):
        """It may setup data structures need by the trajectory. Need not be implemented."""
        self._initialized_read = True

    def write_init(self, system):
        """It may setup data structures need by the trajectory. Need not be implemented."""
        self._initialized_write = True

    def read_sample(self, index): 
        """It must return the sample (system) with the given index"""
        raise NotImplementedError()
        
    def write_sample(self, system, step):
        """It must write a sample (system) to disk. Noting to return."""
        raise NotImplementedError()

    # To read/write timestep and block period sublcasses must
    # implement these methods. They must set the private variable
    # _timestep and _block_period since these properties are cached.
    def read_timestep(self): 
        raise NotImplementedError()

    def write_timestep(self, value): 
        raise NotImplementedError()

    def read_blockperiod(self): 
        raise NotImplementedError()

    def write_blockperiod(self, value): 
        raise NotImplementedError()

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
        # Block period is cached
        if not self._block_period is None:
            return self._block_period

        if len(self.steps) < 2:
            return 1
        delta_old = 0
        delta_one = self.steps[1] - self.steps[0]
        iold = self.steps[0]
        period = 0
        for ii in range(1, len(self.steps)):
            i = self.steps[ii]
            delta = i-iold
            if delta < delta_old and delta == delta_one and abs(delta-delta_old)>2:
                return period
            else:
                period += 1
                iold = i
                delta_old = delta
        return 1

    @block_period.setter
    def block_period(self, value):
        self._block_period = value
        self.write_blockperiod(value)

    def _check_block_period(self):
        """Perform some consistency checks on peridicity of non linear sampling."""
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
        self._unfold()
        msd_total = numpy.sum((self._pos_unf[-1] - self._pos_unf[0])**2) / self._pos_unf[0].shape[0]
        return min(1.0, msd_target * sigma**2 / msd_total) * self.steps[-1] * self.timestep

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

    def split(self, selection=slice(None), index='step', archive=False):
        """Split the trajectory into independent trajectory files, one per sample."""
        if archive:
            tar = tarfile.open(self.filename + '.tar.gz', "w:gz")
        base, ext = os.path.splitext(self.filename)

        for sample, step in zip(self.samples[selection], self.steps[selection]):
            if index == 'step':
                filename = '%s-%09i%s' % (base, step, ext)
            elif index == 'sample':
                filename = '%s-%09i%s' % (base, sample, ext)
            else:
                raise ValueError('unknown option %s' % index)
            t = self.__class__(filename, mode='w')
            t.write_initial_state(self.read_initial_state())
            t.write_sample(self.read_sample(sample), step, sample)
            t.copy_sample(self, step, sample)
            t.close()
            if archive:
                tar.add(filename)
                os.remove(filename)
                
        if archive:
            tar.close()

    # TODO: which one is the good convert()?

    def _convert(self, input_trj, tag=''):
        """Convert trajectory *input_trj* to a different format."""
        # Design: avoid extension driven conversion, only rely on classes.
        # We might pass a prefix or tag to be prepended to self.suffix
        # as in deprecated convert()

        filename = os.path.splitext(self.filename)[0] + tag + '.' + self.suffix
        conv = fmt(filename, mode='w')
        conv.write_initial_state(self.read_initial_state())
        conv.write_timestep(self.timestep)
        conv.write_blockperiod(self.block_period)
        # Careful here: the sample index for writing is not necessarily the same
        # as the one in the inut file, since we might have pruned some cfgs.
        # The initial state is included as the first sample (sample=0)
        for sample, step in zip(self.samples[traj_slice], self.steps[traj_slice]):
            conv.write_sample(self.read_sample(sample), step, sample, ignore=ignore)
        conv.close()
        
        # Return the trajectory object we just closed. 
        # This way we allow post hooks in close().
        # This may slow down things a bit (with xyz format conversion)
        return fmt(filename)
        
    def convert(self, fmt, traj_slice=slice(None), tag='', ignore=[]):

        """
        Convert trajectory object to a different format. 
        *tag* is a prefix to be added before .xyz
        """

        # TODO: convert should not return the converted object. This is meant just to convert the files. 
        # As long as the trajectories implement the interface, all of them should be interchangeable.

        if len(tag)>0:
            tag = '.' + tag

        if isinstance(fmt, str):
            filename = os.path.splitext(self.filename)[0] + tag + '.' + fmt
            if filename == self.filename:
                raise ValueError('Trying to convert a file into itself')
            conv = Trajectory(filename, mode='w', fmt=fmt)
        else:
            # This must be some subclass of Trajectory
            filename = os.path.splitext(self.filename)[0] + tag + '.' + fmt.suffix
            conv = fmt(filename, mode='w')
        conv.write_initial_state(self.read_initial_state())
        # If file has no samples we write the initial state as a sample
        # this should be improved when initial state will be abandoned
        # in favor of some write metadata + write sample.
        if len(self.samples) == 0:
            conv.write_sample(self.read_initial_state(), 0, 0, ignore=ignore)
            
        conv.write_timestep(self.timestep)
        conv.write_blockperiod(self.block_period)
        # Careful here: the sample index for writing is not necessarily the same
        # as the one in the inut file, since we might have pruned some cfgs.
        # The initial state is included as the first sample (sample=0)
        for sample, step in zip(self.samples[traj_slice], self.steps[traj_slice]):
            s = self.read_sample(sample)
            conv.write_sample(s, step, sample, ignore=ignore)
        conv.close()
        
        # Return the trajectory object we just closed. 
        # This way we allow post hooks in close().
        # This may slow down things a bit (with xyz format conversion)
        if isinstance(fmt, str):
            return Trajectory(filename, fmt=fmt)
        else:
            return fmt(filename)


class SuperTrajectory(TrajectoryBase):

    """ Collection of subtrajectories """

    # The approach is inefficient for large number of
    # subtrajectories (init takes some time)

    def __init__(self, subtrajectories, mode='r', timestep=1.0, variable=False):
        f = os.path.dirname(subtrajectories[0].filename) + '/'
        super(SuperTrajectory, self,).__init__(f, mode)
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

        # Redefine samples 
        self.samples = range(len(self.steps))

        # Now it's better check the block period
        self._check_block_period()

    def read_initial_state(self):
        i, j = self._map[0]
        return self.subtrajectories[i].read_sample(j)
                     
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
        
    return __factory_map.get(suffix)(filename, action)
#    return __factory_map.get(suffix)(filename, **kwargs)
