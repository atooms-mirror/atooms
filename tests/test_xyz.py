#!/usr/bin/env python

import unittest
import numpy

from atooms.core.utils import mkdir
from atooms.trajectory import Unfolded
from atooms.trajectory import TrajectoryXYZ, TrajectorySimpleXYZ, TrajectoryRUMD


class TestXYZ(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):
        mkdir('/tmp/test_xyz')
        self.finp = '/tmp/test_xyz/pbc.xyz'
        with open(self.finp, 'w') as fh:
            fh.write("""\
2
1
A 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
2
A 1.1 -1.1 0.0
A -2.9 -2.9 0.0
2
3
A 1.2 -1.2 0.0
A -2.9 2.9 0.0
2
4
A 1.3 -1.3 0.0
A -2.8 2.8 0.0
6.0 6.0 6.0
""")
        # Test metadata recognition
        self.finp_meta = '/tmp/test_xyz/meta.xyz'
        with open(self.finp_meta, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0 step:1 cell:5.0,5.0,5.0 
A 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        # Test 4-dim file
        self.finp_4d = '/tmp/test_xyz/4d.xyz'
        with open(self.finp_4d, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,pos ndim:4 mass:1.0,2.0 step:1
A 1.0 -1.0 0.0 1.0
B 2.9 -2.9 0.0 2.0
2
metafmt:space,comma columns:id,pos mass:1.0,2.0 step:1 cell:5.0,5.0,5.0,5.0
A 1.0 -1.0 0.0 1.0
B 2.9 -2.9 0.0 2.0
2
metafmt:space,comma columns:id,pos mass:1.0,2.0 step:1 cell:5.0,5.0,5.0
A 1.0 -1.0 0.0 1.0
B 2.9 -2.9 0.0 2.0
""")

    def test_xyz_4d(self):
        with self.Trajectory(self.finp_4d) as t:
            meta = t._read_metadata(0)
            self.assertEqual(meta['ndim'], 4)
            self.assertEqual(len(t[0].particle[0].position), 4)
            self.assertEqual(t[0].particle[1].position[3], 2.0)
            # Test to check we grab ndim from the n. of cell coordinates
            meta = t._read_metadata(1)
            self.assertEqual(meta['ndim'], 4)
            self.assertEqual(len(t[1].particle[0].position), 4)
            self.assertEqual(t[1].particle[1].position[3], 2.0)
            # This last test shows that when cell has 3 coords, we are back to ndim=3
            meta = t._read_metadata(2)
            self.assertEqual(meta['ndim'], 3)
            self.assertEqual(len(t[2].particle[0].position), 3)

    def test_xyz_meta(self):
        with self.Trajectory(self.finp_meta) as t:
            meta = t._read_metadata(0)
            self.assertEqual(t.steps, [1])
            self.assertEqual(meta['mass'], [1, 2])
            self.assertEqual(t[0].particle[0].mass, 1.0)
            self.assertEqual(t[0].particle[1].mass, 2.0)

    def test_xyz_indexed(self):
        r_ref = [[1., -1., 0.], [2.9, -2.9, 0.]]
        with self.Trajectory(self.finp) as t:
            self.assertEqual(t.steps, [1, 2, 3, 4])
            self.assertEqual(r_ref[0], list(t[0].particle[0].position))
            self.assertEqual(r_ref[1], list(t[0].particle[1].position))

    def test_xyz_indexed_unfolded(self):
        t1 = self.Trajectory(self.finp)
        t = Unfolded(t1)
        t[0]
        t[1]
        self.assertEqual(list(t[2].particle[1].position), [3.1, -3.1, 0.0])
        self.assertEqual(list(t1[2].particle[1].position), [-2.9, 2.9, 0.0])
        t1.close()
        t.close()

    def test_xyz_indexed_unfolded_skip(self):
        """Test that unfolded trajectories can skip samples"""
        # Here we read all samples
        t1 = self.Trajectory(self.finp)
        t = Unfolded(t1)
        t[0]
        t[1]
        t[2]
        s = t[3]
        t1.close()
        t.close()

        # Here we skip one sample
        t1 = self.Trajectory(self.finp)
        t1 = Unfolded(t1)
        t1[0]
        s1 = t1[3]
        t1.close()

        self.assertEqual(list(s.particle[0].position), list(s1.particle[0].position))
        self.assertEqual(list(s.particle[1].position), list(s1.particle[1].position))

    def test_xyz_unfolded(self):
        """Test reading unfolded after reading normal"""
        with self.Trajectory(self.finp) as t1:
            for s in t1:
                pass
            with Unfolded(t1) as t:
                for s in t:
                    pass

    def test_xyz_with(self):
        r_ref = [[1., -1., 0.], [2.9, -2.9, 0.]]
        with self.Trajectory(self.finp) as t:
            self.assertEqual(r_ref[0], list(t[0].particle[0].position))
            self.assertEqual(r_ref[1], list(t[0].particle[1].position))

        with self.Trajectory(self.finp + '.out', 'w') as to:
            with self.Trajectory(self.finp) as t:
                # This is necessary because timestep is not copied automatically
                # and if the class tries to read it, it won't find
                to.timestep = t.timestep
                to.write(t[0], 0)

        with self.Trajectory(self.finp + '.out') as to:
            self.assertEqual(r_ref[0], list(to[0].particle[0].position))

    def test_xyz_iter(self):
        with self.Trajectory(self.finp) as t:
            for s in t:
                s.particle
            t[-1]

    def test_xyz_to_ram(self):
        # TODO: move this test somewhere else
        from atooms.trajectory import TrajectoryRam
        with self.Trajectory(self.finp) as t:
            with TrajectoryRam() as tram:
                for i, system in enumerate(t):
                    tram[i] = system
                    # Write a second time, this will overwrite
                    tram[i] = system
                    self.assertEqual(tram[i].particle[0].position[0],
                                     system.particle[0].position[0])

    def test_xyz_ids(self):
        # Test ids ordering
        finp = '/tmp/test_xyz/meta.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z step:1 cell:5.0,5.0,5.0 
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z step:2 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            self.assertEqual(th[0].particle[0].species, 'B')
            self.assertEqual(th[0].particle[1].species, 'A')
            self.assertEqual(th[1].particle[0].species, 'C')
            self.assertEqual(th[1].particle[1].species, 'B')

    def test_xyz_mass(self):
        finp = '/tmp/test_xyz/meta.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
3
metafmt:space,comma columns:id,x,y,z,radius mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0 
B 1.0 -1.0 0.0 0.5
A 2.9 -2.9 0.0 0.6
C 2.9 -2.9 0.0 0.7
3
metafmt:space,comma columns:id,x,y,z mass:2.0,3.0 step:1 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
B 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            self.assertEqual(th[0].particle[0].radius, 0.5)
            self.assertEqual(th[0].particle[1].radius, 0.6)
            self.assertEqual(th[0].particle[2].radius, 0.7)
            self.assertEqual(th[0].particle[0].mass, 2.0)
            self.assertEqual(th[0].particle[1].mass, 1.0)
            self.assertEqual(th[0].particle[2].mass, 3.0)
            self.assertEqual(th[1].particle[0].mass, 3.0)
            self.assertEqual(th[1].particle[1].mass, 2.0)
            self.assertEqual(th[1].particle[2].mass, 2.0)

        # This is ok
        with open(finp, 'w') as fh:
            fh.write("""\
3
columns:id,x,y,z mass:2 step:1
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
C 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            self.assertEqual(th[0].particle[0].mass, 2.0)
            self.assertEqual(th[0].particle[1].mass, 2.0)
            self.assertEqual(th[0].particle[2].mass, 2.0)

        # This is not ok, mass / species mismatch
        with open(finp, 'w') as fh:
            fh.write("""\
3
columns:id,x,y,z mass:1,2 step:1
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
C 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            with self.assertRaises(ValueError):
                th[0]

    def test_xyz_columns(self):
        finp = '/tmp/test_xyz/columns.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
1
columns:name,pos,vel step:1
A 1.0 -1.0 0.0 1.0 2.0 3.0
""")
        with self.Trajectory(finp) as th:
            self.assertEqual(th._read_metadata(0)['columns'], ['name', 'pos', 'vel'])
            self.assertEqual(list(th[0].particle[0].velocity), [1.0, 2.0, 3.0])

    def test_xyz_cell(self):
        """Test that fluctuating cell side is read correctly."""
        finp = '/tmp/test_xyz/cell.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
1
step:0 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:1 cell:2.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:2 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:3 cell:3.0,1.0,1.0
A 1.0 -1.0 0.0
""")
        with TrajectoryXYZ(finp) as th:
            self.assertEqual([s.cell.side[0] for s in th], [1.0, 2.0, 1.0, 3.0])

    def tearDown(self):
        from atooms.core.utils import rmd
        rmd('/tmp/test_xyz')


from atooms.trajectory.xyz import TrajectoryNeighbors

class TestNeighbors(unittest.TestCase):

    Trajectory = TrajectoryNeighbors

    def setUp(self):
        mkdir('/tmp/test_xyz')

    def test_neighbors_consume(self):
        with open('/tmp/test_xyz/neighbors.xyz', 'w') as fh:
            fh.write("""\
4
step:1 columns:neighbors*
2 4
1 3
2 
1
""")
        with TrajectoryNeighbors('/tmp/test_xyz/neighbors.xyz', offset=0) as th:
            s = th[0]
            self.assertEqual(list(s.particle[0].neighbors), [2, 4])
            self.assertEqual(list(s.particle[1].neighbors), [1, 3])
            self.assertEqual(list(s.particle[2].neighbors), [2])
            self.assertEqual(list(s.particle[3].neighbors), [1])

    def test_neighbors_comma(self):
        with open('/tmp/test_xyz/neighbors.xyz', 'w') as fh:
            fh.write("""\
4
step:1 columns:neighbors timestep:0.001
2,4
1,3
2
1
""")
        with TrajectoryNeighbors('/tmp/test_xyz/neighbors.xyz', offset=0) as th:
            s = th[0]
            self.assertEqual(th.timestep, 0.001)
            self.assertEqual(list(s.particle[0].neighbors), [2, 4])
            self.assertEqual(list(s.particle[1].neighbors), [1, 3])
            self.assertEqual(list(s.particle[2].neighbors), [2])
            self.assertEqual(list(s.particle[3].neighbors), [1])

    def tearDown(self):
        from atooms.core.utils import rmd
        rmd('/tmp/test_xyz')


class TestSimpleXYZ(TestXYZ):

    Trajectory = TrajectorySimpleXYZ

    def test_xyz_4d(self):
        pass

    def test_xyz_meta(self):
        pass

    def test_xyz_mass(self):
        pass

    def test_xyz_columns(self):
        pass

    def test_xyz_cell(self):
        pass


class TestRumd(unittest.TestCase):

    Trajectory = TrajectoryRUMD

    def setUp(self):
        mkdir('/tmp/test_xyz')
        super(TestRumd, self).setUp()
        ioformat_1 = """\
2
ioformat=1 dt=0.005000000 boxLengths=6.000000000,6.000000000,6.000000000 numTypes=2 Nose-Hoover-Ps=-0.154678583 Barostat-Pv=0.000000000 mass=1.000000000,1.000000000 columns=type,x,y,z,imx,imy,imz,vx,vy,vz,fx,fy,fz,pe,vir
0  1.0 -1.0 0.0 0 0 1 -0.955896854 -2.176721811 0.771060944 14.875996590 -28.476327896 -15.786120415 -5.331668218 22.538120270
1  2.9 -2.9 0.0 0 0 -1 -0.717318892 -0.734904408 0.904034972 -28.532371521 13.714955330 0.387423307 -7.276362737 11.813765526
"""

        with open('/tmp/test_xyz/rumd.xyz', 'w') as fh:
            fh.write(ioformat_1)
            self.finp = fh.name

    def test_read_write(self):
        with TrajectoryRUMD(self.finp) as th:
            with TrajectoryRUMD(self.finp + 'out', 'w') as th_out:
                self.assertEqual(th.timestep, 0.005)
                self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
                th_out.timestep = th.timestep
                for i, system in enumerate(th):
                    th_out.write(system, th.steps[i])
        with TrajectoryRUMD(self.finp + 'out') as th:
            self.assertEqual(th.timestep, 0.005)
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])


class TestUtils(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):
        mkdir('/tmp/test_xyz')
        with open('/tmp/test_xyz/1.xyz', 'w') as fh:
            fh.write("""\
1
metafmt:space,comma columns:id,x,y,z step:1 cell:5.0,5.0,5.0
B 1.0 -1.0 0.0
1
metafmt:space,comma columns:id,x,y,z step:2 cell:5.0,5.0,5.0 
B 2.9 -2.9 0.0
1
metafmt:space,comma columns:id,x,y,z step:4 cell:5.0,5.0,5.0 
B 2.9 -2.9 0.0
""")
        with open('/tmp/test_xyz/2.xyz', 'w') as fh:
            fh.write("""\
1
metafmt:space,comma columns:id,x,y,z step:2 cell:5.0,5.0,5.0
B 1.0 -1.0 0.0
1
metafmt:space,comma columns:id,x,y,z step:3 cell:5.0,5.0,5.0 
B 2.9 -2.9 0.0
1
metafmt:space,comma columns:id,x,y,z step:4 cell:5.0,5.0,5.0 
B 2.9 -2.9 0.0
1
metafmt:space,comma columns:id,x,y,z step:5 cell:5.0,5.0,5.0 
B 2.9 -2.9 0.0
""")

    def test_paste(self):
        import atooms.trajectory as trj
        t1 = trj.TrajectoryXYZ('/tmp/test_xyz/1.xyz')
        t2 = trj.TrajectoryXYZ('/tmp/test_xyz/2.xyz')
        steps = []
        for step, s1, s2 in trj.utils.paste(t1, t2):
            steps.append(step)
        self.assertEqual(steps, [2, 4])
        t1.close()
        t2.close()

    def tearDown(self):
        from atooms.core.utils import rmd
        rmd('/tmp/test_xyz')


if __name__ == '__main__':
    unittest.main()
