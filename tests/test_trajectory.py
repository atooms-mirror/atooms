#!/usr/bin/env python

import os
import unittest
import numpy

from atooms.core.utils import rmd, rmf
from atooms.system import System, Particle, Cell
from atooms.trajectory import Unfolded
from atooms.trajectory import TrajectoryXYZ, TrajectorySimpleXYZ, TrajectoryRUMD
import atooms.trajectory as trj


def _equal(system1, system2, ignore=None):
    check = {}
    check['npart'] = len(system1.particle) == len(system2.particle)
    for p1, p2 in zip(system1.particle, system2.particle):
        check['position'] = all(p1.position == p2.position)
        check['mass'] = p1.mass == p2.mass
        check['species'] = p1.species == p2.species
    check['side'] = all(system1.cell.side == system1.cell.side)
    for key in check:
        if ignore is not None and key in ignore:
            continue
        if not check[key]:
            print key, 'differs'
            return False
    return True


def _rename_species(particle, db):
    for p in particle:
        p.species = db[p.species]
    return particle


class Test(unittest.TestCase):

    def setUp(self):
        import copy
        particle = [Particle(position=[0.0, 0.0, 0.0], species='A', mass=1.0),
                    Particle(position=[1.0, 1.0, 1.0], species='B', mass=2.0),
                    ]
        cell = Cell([2.0, 2.0, 2.0])
        self.system = []
        self.system.append(System(copy.deepcopy(particle), cell))
        self.system.append(System(copy.deepcopy(particle), cell))
        self.inpfile = '/tmp/test_trajectory'
        self.inpdir = '/tmp/test_trajectory.d'
        from atooms.core.utils import mkdir
        mkdir(self.inpdir)

    def _read_write(self, cls, path=None, ignore=None):
        """Read and write"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            th.write_timestep(1.0)
            for i, system in enumerate(self.system):
                th.write(self.system[i], i)
        with cls(path) as th:
            self.assertEqual(th.timestep, 1.0)
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore))

    def _write(self, cls, path=None):
        """Write only"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            th.write_timestep(1.0)
            for i, system in enumerate(self.system):
                th.write(self.system[i], i)

    def test_xyz(self):
        # TODO: mass is not written by xyz
        self._read_write(trj.TrajectoryXYZ, ignore=['mass'])
        self._read_write(trj.TrajectorySimpleXYZ, ignore=['mass'])

    def test_ram(self):
        self._read_write(trj.TrajectoryRam)
        self._read_write(trj.ram.TrajectoryRamFull)

    def test_hdf5(self):
        self._read_write(trj.TrajectoryHDF5)

    def test_rumd(self):
        # RUMD uses integer ids for checmical species. They should be integers.
        for s in self.system:
            s.particle = _rename_species(s.particle, {'A': '0', 'B': '1'})
        self._read_write(trj.TrajectoryRUMD)
        # TODO: add write_sample() to supertrajectory
        #self._read_write(trj.SuperTrajectoryRUMD, self.inpdir, ignore=['id', 'name'])

    def test_pdb(self):
        self._write(trj.TrajectoryPDB)
        reference = """\
MODEL        0
CRYST1    2.000    2.000    2.000     90     90     90 P 1           1
HETATM    0             A       0.000   0.000   0.000  1.00  1.00             A
HETATM    1             B       1.000   1.000   1.000  1.00  1.00             B
MODEL        1
CRYST1    2.000    2.000    2.000     90     90     90 P 1           1
HETATM    0             A       0.000   0.000   0.000  1.00  1.00             A
HETATM    1             B       1.000   1.000   1.000  1.00  1.00             B
"""
        output = open(self.inpfile).read()
        self.assertTrue(output == reference)

    def tearDown(self):
        rmf(self.inpfile)
        rmd(self.inpdir)


if __name__ == '__main__':
    unittest.main()
