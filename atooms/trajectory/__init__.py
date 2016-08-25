# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

import os
import sys
import numpy
import tarfile
import logging
import warnings

from decorators import *
from base import TrajectoryBase

# Useful functions to manipulate trajectories

def convert(inp, out, fout='', tag='', prefix='', force=True, exclude=[], include=[], stdout=False, callback=None, args={}):
    """Convert trajectory into a different format.

    inp: input trajectory object
    out: output trajectory class
    tag: optional string to be prepended before the output suffix

    Return: name of converted trajectory file
    """
    # TODO: convert metadata (interaction etc) !
    # If the input trajectory lies in a directory, the new trajectory is located
    # in a companion directory prefixed by tag. The basename is config

    # Check that we have some files there
    if len(inp) == 0:
        raise IOError('no samples in trajectory (%s)' % inp.filename)

    if stdout:
        filename = '/dev/stdout'
    else:
        if len(fout) > 0:
            filename = fout
        elif os.path.isdir(inp.filename):
            if tag == '':
                if prefix == '':
                    tag = '-conv'
            d = os.path.dirname(inp.filename)
            b = os.path.basename(inp.filename)
            dirname = os.path.join(d, prefix + b) + tag
            from pyutils.utils import mkdir
            mkdir(dirname)
            filename = dirname + '/config.' + out.suffix
        else:
            filename = os.path.splitext(inp.filename)[0] + tag + '.' + out.suffix    

    if not os.path.exists(filename) or force:
        with out(filename, 'w') as conv:
            conv.exclude(exclude)
            conv.include(include)
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

from base import TrajectoryBase

class SuperTrajectory(TrajectoryBase):

    """ Collection of subtrajectories """

    # Ppproach inefficient for large number of files (init takes some time)
    # The path of a supertrajectory should be the directory containing all the files? Mmmhh

    def __init__(self, files, trajectoryclass, mode='r', variable=False, periodic=True):
        self.subtrajectories = [trajectoryclass(f, mode=mode) for f in files]
        f = os.path.dirname(self.subtrajectories[0].filename) + '/'
        super(SuperTrajectory, self).__init__(f, mode)
        self._timestep = 1.0
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

            # TODO: drop zipping of steps and samples
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


class SuperTrajectory2(TrajectoryBase):

    """Collection of subtrajectories"""

    # Optimized version

    def __init__(self, files, trajectoryclass, mode='r'):
        """Group list of files into a single trajectory"""
        self.files = files
        if len(self.files) == 0:
            raise ValueError('no files found in %s' % self.files)
        f = os.path.dirname(self.files[0])
        super(SuperTrajectory2, self).__init__(f, mode)
        self.trajectoryclass = trajectoryclass
        self._timestep = 1.0

        # Make sure subtrajectories are sorted by increasing step
        self.files.sort()
        self._map = []
        for i, f in enumerate(self.files):
            # This is slow... just to get the step index....
            with self.trajectoryclass(f) as t:
                self._map.append((t.steps[0], f))
            #self._map.append((i, f))
        self.steps = [i[0] for i in self._map]

    def read_sample(self, sample):
        step, f = self._map[sample]
        with self.trajectoryclass(f) as t:
            return t[0]


# Factory method which mimics an abstract factory class
__factory_map = {}

from xyz import *
from pdb import TrajectoryPDB
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

# Load plugins (if plugins are found)
try:
    from atooms.plugins.trajectory import *
except:
    # No plugins found
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
