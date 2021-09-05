#!/usr/bin/env python

import os
import unittest
import numpy

from atooms.core.utils import rmd, rmf
from atooms.system import System, Particle, Cell
from atooms.trajectory import Unfolded
from atooms.trajectory import TrajectoryXYZ, TrajectorySimpleXYZ, TrajectoryRUMD, Trajectory
from atooms.trajectory.base import TrajectoryBase
import atooms.trajectory as trj

def almost_equal(x, y, rtol):
    try:
        numpy.testing.assert_allclose(x, y, rtol=rtol)
        return True
    except AssertionError:
        return False

def _equal(system1, system2, ignore=None, verbose=True, precision=1e-10):
    check = {}
    check['npart'] = len(system1.particle) == len(system2.particle)
    check['side'] = all(system1.cell.side == system1.cell.side)
    for p1, p2 in zip(system1.particle, system2.particle):
        # There is roundoffs with gsd! It is surprising because it is a binary format. It should be up to machine precision... difference is 10 places
        # check['position'] = all(p1.position == p2.position)
        check['position'] = almost_equal(p1.position, p2.position, rtol=precision)
        check['velocity'] = almost_equal(p1.velocity, p2.velocity, rtol=precision)
        check['mass'] = p1.mass == p2.mass
        check['species'] = p1.species == p2.species
    for key in check:
        if ignore is not None and key in ignore:
            continue
        if not check[key]:
            if verbose:
                print(key, 'differs')
            return False
    return True

def _difference(system1, system2, ignore=None, precision=1e-10):
    check = {}
    check['npart'] = len(system1.particle) == len(system2.particle)
    check['side'] = all(system1.cell.side == system1.cell.side)
    for p1, p2 in zip(system1.particle, system2.particle):
        check['position'] = almost_equal(p1.position, p2.position, rtol=precision)
        check['velocity'] = almost_equal(p1.velocity, p2.velocity, rtol=precision)
        check['mass'] = p1.mass == p2.mass
        check['species'] = p1.species == p2.species
    diffs = []
    for key in check:
        if ignore is not None and key in ignore:
            continue
        if not check[key]:
            diffs.append(key)
    return diffs


def _rename_species(particle, db):
    for p in particle:
        p.species = db[p.species]
    return particle


class Test(unittest.TestCase):

    def setUp(self):
        import copy
        particle = [Particle(position=[0.1, 0.2, 0.3], velocity=[0.3, 0.2, 0.1], species='A', mass=1.0),
                    Particle(position=[1.3, 1.2, 1.1], velocity=[1.1, 1.2, 1.3], species='B', mass=2.0),
                    ]
        cell = Cell([2.0, 2.0, 2.0])
        self.system = []
        self.system.append(System(copy.deepcopy(particle), cell))
        self.system.append(System(copy.deepcopy(particle), cell))
        self.inpfile = '/tmp/test_trajectory'
        self.inpdir = '/tmp/test_trajectory.d'
        from atooms.core.utils import mkdir
        mkdir(self.inpdir)

    def _read_write(self, cls, path=None, ignore=None, precision=1e-10):
        """Read and write"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            th.write_timestep(1.0)
            for i, system in enumerate(self.system):
                th.write(system, i)
        with cls(path) as th:
            self.assertEqual(th.timestep, 1.0)
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore, precision=precision))
                self.assertTrue(self.system[i].__class__ is system.__class__)

    def _read_write_fields(self, cls, write_fields=None, read_fields=None, path=None, ignore=None, fail=None, precision=1e-10):
        """Read and write with fields"""
        if path is None:
            path = self.inpfile

        # Write
        # try:
        #     th = cls(path, 'w', fields=write_fields)
        # except TypeError:
        th = cls(path, 'w')
        th.variables = write_fields
        th.write_timestep(1.0)
        for i, system in enumerate(self.system):
            th.write(system, i)
        th.close()

        # Read
        # try:
        #     th = cls(path, fields=read_fields)
        # except TypeError:
        th = cls(path)
        th.variables = read_fields
        self.assertEqual(th.timestep, 1.0)
        for i, system in enumerate(th):
            if fail is not None:
                self.assertTrue(set(_difference(self.system[i], system, ignore, precision=precision)) == set(fail))
            else:
                self.assertTrue(_equal(self.system[i], system, ignore, verbose=False, precision=precision))
            self.assertTrue(self.system[i].__class__ is system.__class__)
        th.close()

    def _write(self, cls, path=None):
        """Write only"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            th.write_timestep(1.0)
            for i, system in enumerate(self.system):
                th.write(self.system[i], i)

    def _slice(self, cls, path=None):
        """Write only"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            for i, system in enumerate(self.system):
                th.write(system, i)
            from atooms.trajectory.decorators import Sliced
            ts = Sliced(th, slice(None, None, 2))
            self.assertEqual(len(th), 2)
            self.assertEqual(len(ts), 1)
            # This will fail, we cannot slice twice yet
            # ts = Sliced(ts, slice(None, None, 1))
            # print len(ts)

    def _append(self, cls, path=None, ignore=None):
        """Read and write"""
        if path is None:
            path = self.inpfile
        with cls(path, 'w') as th:
            for system in self.system:
                th.append(system)
        with cls(path) as th:
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore))
                self.assertTrue(self.system[i].__class__ is system.__class__)

    def _convert(self, cls_inp, cls_out, path=None, ignore=None):
        """Write then convert"""
        if path is None:
            path = self.inpfile
        fout = self.inpfile + '.out'
        with cls_inp(path, 'w') as th:
            th.write_timestep(1.0)
            for i, system in enumerate(self.system):
                th.write(self.system[i], i)

        with cls_inp(path) as th:
            from atooms.trajectory.utils import convert
            _ = convert(th, cls_out, fout)

        if isinstance(cls_out, str):
            th = trj.Trajectory(fout, fmt=cls_out)
        else:
            th = cls_out(fout)
        self.assertEqual(th.timestep, 1.0)
        self.assertTrue(len(th.steps), len(self.system))
        for i, system in enumerate(th):
            self.assertTrue(_equal(self.system[i], system, ignore))
        th.close()

    def test_xyz(self):
        # TODO: mass is not written by xyz
        # Velocity is not written in default xyz format, so we ignore it
        self._read_write(trj.TrajectoryXYZ, ignore=['mass', 'velocity'])
        self._convert(trj.TrajectoryXYZ, trj.TrajectoryXYZ, ignore=['mass', 'velocity'])
        self._convert(trj.TrajectoryXYZ, 'xyz', ignore=['mass', 'velocity'])
        self._append(trj.TrajectoryXYZ, ignore=['mass', 'velocity'])
        self._slice(trj.TrajectoryXYZ)
        # Check that when requesting to write the velocity, we actually write it and read it back automatically
        self._read_write_fields(trj.TrajectoryXYZ, write_fields=['species', 'position', 'velocity'], ignore=['mass'])
        # This must fail: writing positions but not reading them
        self._read_write_fields(trj.TrajectoryXYZ, write_fields=['species', 'position', 'velocity'], read_fields=['species', 'velocity'], ignore=['mass'], fail=['position'])

    def test_simple_xyz(self):
        # Mass and velocity is not written in default xyz format, so we ignore it
        self._read_write(trj.TrajectorySimpleXYZ, ignore=['mass', 'velocity'])
        self._append(trj.TrajectorySimpleXYZ, ignore=['mass', 'velocity'])

    def test_ram(self):
        self._read_write(trj.TrajectoryRam)
        self._append(trj.TrajectoryRam)

    def test_ram_full(self):
        self._read_write(trj.ram.TrajectoryRamFull)
        self._append(trj.ram.TrajectoryRamFull)

    def test_hdf5(self):
        try:
            import h5py
        except ImportError:
            self.skipTest('missing hdf5')

        self._read_write(trj.TrajectoryHDF5)
        self._read_write_fields(trj.TrajectoryHDF5, write_fields=['species', 'position', 'velocity'], read_fields=['species', 'position', 'velocity'])
        # This must fail: writing velocities but not reading them
        self._read_write_fields(trj.TrajectoryHDF5, write_fields=['species', 'position'], read_fields=['species', 'position', 'velocity'], fail=['velocity'])
        # Velocity is not kept in conversion to xyz
        self._convert(trj.TrajectoryXYZ, 'hdf5', ignore=['mass', 'velocity'])

    def test_gsd(self):
        try:
            import gsd
            from atooms.trajectory.gsd import TrajectoryGSD
        except ImportError:
            self.skipTest('missing gsd')

        self._read_write(trj.TrajectoryGSD, ignore=['mass', 'velocity'], precision=1e-7)
        #self._convert(trj.TrajectoryGSD, 'gsd', ignore=['mass', 'velocity'])
        #self._read_write_fields(trj.TrajectoryGSD, write_fields=['species', 'position', 'velocity'], read_fields=['species', 'position', 'velocity'], ignore=['mass'])
        ## This must fail: writing velocities but not reading them
        #self._read_write_fields(trj.TrajectoryGSD, write_fields=['species'], read_fields=['species', 'position'], ignore=['mass'], fail=['velocity'])

    def test_rumd(self):
        # RUMD uses integer ids for chemical species. They should be integers.
        for s in self.system:
            s.particle = _rename_species(s.particle, {'A': '0', 'B': '1'})
        self._read_write(trj.TrajectoryRUMD)
        # TODO: all these fail without changing _optimize_fields() in setup_fields() so as to return a new list
        # self._read_write_fields(trj.TrajectoryRUMD, read_fields=['species', 'position', 'velocity'], ignore=['mass'])
        # self._read_write_fields(trj.TrajectoryRUMD, read_fields=['type', 'x', 'y', 'z', 'vx', 'vy', 'vz'], ignore=['mass'])
        # This must fail: writing velocities but not reading them
        # self._read_write_fields(trj.TrajectoryRUMD, read_fields=['species', 'position'], ignore=['mass'], fail=True)
        # self._read_write_fields(trj.TrajectoryRUMD, read_fields=['type', 'x', 'y', 'z'], ignore=['mass'], fail=True)
        # TODO: add write_sample() to supertrajectory
        #self._read_write(trj.SuperTrajectoryRUMD, self.inpdir, ignore=['id', 'name'])

    def test_pdb(self):
        self._read_write(trj.TrajectoryPDB, ignore=['mass', 'velocity'])
        with trj.TrajectoryPDB('data/trajectory.pdb') as th:
            self.assertTrue(numpy.all(th[0].cell.side == numpy.array([10.0, 10.0, 10.0])))
            self.assertTrue(numpy.all(th[1].cell.side == numpy.array([10.0, 10.0, 10.0])))
        
    def test_super(self):
        import glob
        with TrajectoryXYZ(os.path.join(self.inpdir, '0.xyz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 10)
        with TrajectoryXYZ(os.path.join(self.inpdir, '1.xyz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 20)
        with trj.SuperTrajectory(glob.glob(self.inpdir + '/*'), TrajectoryXYZ) as th:
            self.assertFalse(th.grandcanonical)
            self.assertEqual(th.times, [0.001*10, 0.001*20])
            self.assertEqual(th.timestep, 0.001)
            self.assertEqual(th.steps, [10, 20])
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore=['mass', 'velocity']))

    def test_super_rumd(self):
        # TODO: refactor supertrajectory tests
        # RUMD uses integer ids for chemical species. They should be integers.
        for s in self.system:
            s.particle = _rename_species(s.particle, {'A': '0', 'B': '1'})
        with trj.TrajectoryRUMD(os.path.join(self.inpdir, '0.xyz.gz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 0)
        with trj.TrajectoryRUMD(os.path.join(self.inpdir, '1.xyz.gz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 1)
        with trj.rumd.SuperTrajectoryRUMD(self.inpdir) as th:
            self.assertEqual(th.times, [0.001*0, 0.001*1])
            self.assertEqual(th.timestep, 0.001)
            self.assertEqual(th.steps, [0, 1])
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore=['mass']))

    @unittest.expectedFailure
    def test_super_rumd_fields(self):
        for s in self.system:
            s.particle = _rename_species(s.particle, {'A': '0', 'B': '1'})
        with trj.TrajectoryRUMD(os.path.join(self.inpdir, '0.xyz.gz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 0)
        with trj.TrajectoryRUMD(os.path.join(self.inpdir, '1.xyz.gz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 1)
        # TODO: supertrajectories should propagate fields but currently do not
        with trj.rumd.SuperTrajectoryRUMD(self.inpdir) as th:
            th.fields = ['type', 'x', 'y', 'z']
            system = th[0]
            self.assertTrue(_equal(self.system[0], system, ignore=['mass']), fail=True)
            self.assertTrue(self.system[0].__class__ is system.__class__)

    def test_folder(self):
        import glob
        with TrajectoryXYZ(os.path.join(self.inpdir, '10.xyz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 10)
        with TrajectoryXYZ(os.path.join(self.inpdir, '20.xyz'), 'w') as th:
            th.timestep = 0.001
            th.write(self.system[0], 20)
        with trj.folder.Foldered(self.inpdir, cls='xyz') as th:
            self.assertEqual(th.times, [0.001*10, 0.001*20])
            self.assertEqual(th.timestep, 0.001)
            self.assertEqual(th.steps, [10, 20])
            for i, system in enumerate(th):
                self.assertTrue(_equal(self.system[i], system, ignore=['mass', 'velocity']))

    def test_cache(self):
        with TrajectoryXYZ(os.path.join(self.inpdir, 'cache.xyz'), 'w') as th:
            for i in range(10):
                th.write(self.system[0], i)

        with TrajectoryXYZ(os.path.join(self.inpdir, 'cache.xyz'), 'r') as th_0:
            s_0 = th_0[0]

        with TrajectoryXYZ(os.path.join(self.inpdir, 'cache.xyz'), 'r') as th_1:
            th_1.cache = True
            s_1 = th_1.read(0)
            s_1 = th_1.read(0)

    def test_block_size(self):
        from atooms.trajectory.utils import check_block_size, get_block_size
        steps = [2**i for i in range(5)]
        block = get_block_size(steps)
        check_block_size(steps, block)

    def test_decorator(self):
        """Test that frame is accessible to callbacks"""
        finp = os.path.join(self.inpdir, 'test.xyz')
        with open(finp, 'w') as fh:
            fh.write("""\
2
columns:id,x,y,z step:1 cell:5.0,5.0,5.0
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
columns:id,x,y,z step:2 cell:5.0,5.0,5.0
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        def cbk(system):
            system.frame
            return system
        with TrajectoryXYZ(finp) as th:
            th.add_callback(cbk)
            th[0]
            th[1]

    def test_class_callback(self):
        def f(s):
            s._signal = True
            return s

        class TrajectoryXYZCustom(TrajectoryXYZ):
            pass
        TrajectoryXYZCustom.add_class_callback(f)
        Trajectory.add(TrajectoryXYZCustom)

        self.assertFalse(TrajectoryXYZCustom.class_callbacks is TrajectoryXYZ.class_callbacks)

        with Trajectory(os.path.join(self.inpdir, 'test.xyz'), 'w') as th:
            th.write(self.system[0])
        with Trajectory(os.path.join(self.inpdir, 'test.xyz'), 'r') as th:
            s = th[0]
            self.assertTrue(hasattr(s, '_signal'))

        TrajectoryXYZCustom.class_callbacks.remove((f, (), {}))
        with Trajectory(os.path.join(self.inpdir, 'test.xyz'), 'r') as th:
            s = th[0]
            self.assertFalse(hasattr(s, '_signal'))

    def _copy_inplace(self, trajectory, expect=False):
        """
        Test that trajectory returns a copy of the system and that
        modifications are not propagated to the underlying trajectory.

        Test in-place modification
        """
        system = trajectory[0]
        original = system.particle[1].position.copy()
        system.particle[1].position *= 2
        new_system = trajectory[0]
        if expect:
            self.assertEqual(system.particle[1].position[0], new_system.particle[1].position[0])
        else:
            self.assertNotEqual(system.particle[1].position[0], new_system.particle[1].position[0])

    def _copy_reassign(self, trajectory, expect=False):
        """
        Test that trajectory returns a copy of the system and that
        modifications are not propagated to the underlying trajectory.

        Test assignement
        """
        system = trajectory[0]
        original = system.particle[1].position.copy()
        system.particle[1].position[0] = 100000000000000.0
        if expect:
            self.assertEqual(system.particle[1].position[0], new_system.particle[1].position[0])
        else:
            self.assertNotEqual(system.particle[1].position[0], new_system.particle[1].position[0])

    def test_copy_ram_view(self):
        with trj.ram.TrajectoryRamView() as th:
            th[0] = self.system[0]
            self._copy_inplace(th, expect=True)

    def test_copy_ram(self):
        with trj.ram.TrajectoryRam() as th:
            th[0] = self.system[0]
            self._copy_inplace(th)

    def test_copy_xyz(self):
        with trj.TrajectoryXYZ(self.inpfile, 'w') as th:
            th.write(self.system[0])
        with trj.TrajectoryXYZ(self.inpfile, 'r') as th:
            self._copy_inplace(th)

    def test_info(self):
        with trj.xyz.TrajectoryXYZ(self.inpfile, 'w') as th:
            th.write(self.system[0])
            th.write(self.system[1])
        with trj.xyz.TrajectoryXYZ(self.inpfile, 'r') as th:
            from atooms.trajectory.utils import info
            info(th)
            keys = [
                'path',
                'format',
                'frames',
                'megabytes',
                'particles',
                'species',
                'composition',
                'cell density',
                'cell side',
                'cell volume',
                'steps',
                'duration',
                'timestep',
                'block size',
                'steps between frames',
                'time between frames',
                'block steps',
                'block',
                'grandcanonical']
            info(th, keys=','.join(keys))

    def test_lammps(self):
        import sys
        with open(self.inpfile, 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-3 3
-3 3
-3 3
ITEM: ATOMS id type xs ys zs
2 1 0.10 0.11 0.12
1 1 0.20 0.21 0.22
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-4 4
-4 4
-4 4
ITEM: ATOMS id type xs ys zs
1 1 0.00 0.01 0.02
2 1 0.50 0.51 0.52
""")
        from atooms.trajectory import TrajectoryLAMMPS
        def scale(pos, side):
            return [(x - 0.5) * L for x, L in zip(pos, side)]
        with TrajectoryLAMMPS(self.inpfile) as th:
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
            self.assertEqual(list(th[0].particle[0].position), scale([0.20, 0.21, 0.22], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[0].particle[1].position), scale([0.10, 0.11, 0.12], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[1].cell.side), [8.0, 8.0, 8.0])
            self.assertEqual(list(th[1].particle[0].position), scale([0.00, 0.01, 0.02], [8.0, 8.0, 8.0]))
            self.assertEqual(list(th[1].particle[1].position), scale([0.50, 0.51, 0.52], [8.0, 8.0, 8.0]))

    def test_lammps_folder(self):
        with open(self.inpdir + '/0.atom', 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-3 3
-3 3
-3 3
ITEM: ATOMS id type xs ys zs
2 1 0.10 0.11 0.12
1 1 0.20 0.21 0.22
""")
        with open(self.inpdir + '/1.atom', 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
20
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-4 4
-4 4
-4 4
ITEM: ATOMS id type xs ys zs
1 1 0.00 0.01 0.02
2 1 0.50 0.51 0.52
""")
        from atooms.trajectory import TrajectoryFolderLAMMPS
        def scale(pos, side):
            return [(x - 0.5) * L for x, L in zip(pos, side)]
        with TrajectoryFolderLAMMPS(self.inpdir) as th:
            self.assertEqual(th.steps, [10, 20])
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
            self.assertEqual(list(th[0].particle[0].position), scale([0.20, 0.21, 0.22], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[0].particle[1].position), scale([0.10, 0.11, 0.12], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[1].cell.side), [8.0, 8.0, 8.0])
            self.assertEqual(list(th[1].particle[0].position), scale([0.00, 0.01, 0.02], [8.0, 8.0, 8.0]))
            self.assertEqual(list(th[1].particle[1].position), scale([0.50, 0.51, 0.52], [8.0, 8.0, 8.0]))

    def tearDown(self):
        rmf(self.inpfile)
        rmd(self.inpdir)


if __name__ == '__main__':
    unittest.main()
