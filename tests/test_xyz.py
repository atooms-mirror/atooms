#!/usr/bin/env python

import unittest
import numpy

from atooms.trajectory import Unfolded
from atooms.trajectory import TrajectoryXYZ, TrajectorySimpleXYZ, TrajectoryRUMD

class TestXYZ(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):        
        self.finp = '/tmp/test_pbc.xyz'
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
        self.finp_meta = '/tmp/test_meta.xyz'
        with open(self.finp_meta, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0 step:1 cell:5.0,5.0,5.0 
A 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        # Test 4-dim file
        self.finp_4d = '/tmp/test_4d.xyz'
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
        t = self.Trajectory(self.finp)
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
        t.close()

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

    def test_xyz_ids(self):
        # Test ids ordering
        finp = '/tmp/test_meta.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0 
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            th._id_min = 1
            self.assertEqual(th[0].particle[0].id, 2)
            self.assertEqual(th[0].particle[1].id, 1)
            self.assertEqual(th[1].particle[0].id, 3)
            self.assertEqual(th[1].particle[1].id, 2)

    def test_xyz_mass(self):
        finp = '/tmp/test_meta.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0 
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        with self.Trajectory(finp) as th:
            th._id_min = 1
            self.assertEqual(th[0].particle[0].mass, 2.0)
            self.assertEqual(th[0].particle[1].mass, 1.0)
            self.assertEqual(th[1].particle[0].mass, 3.0)
            self.assertEqual(th[1].particle[1].mass, 2.0)

    def test_xyz_columns(self):
        finp = '/tmp/test_xyz_columns.xyz'
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
        finp = '/tmp/test_xyz_cell.xyz' 
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


class TestSimpleXYZ(TestXYZ):

    # TODO: refactor generic tests for trajectories like this :-)
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


class TestRumd(TestXYZ):

    # TODO: refactor generic tests for trajectories like this :-)
    Trajectory = TrajectoryRUMD

    def setUp(self):
        super(TestRumd, self).setUp()
        ioformat_1 = """\
2
ioformat=1 dt=0.005000000 boxLengths=6.000000000,6.000000000,6.000000000 numTypes=2 Nose-Hoover-Ps=-0.154678583 Barostat-Pv=0.000000000 mass=1.000000000,1.000000000 columns=type,x,y,z,imx,imy,imz,vx,vy,vz,fx,fy,fz,pe,vir
0 2.545111895 -0.159052730 -2.589233398 0 0 1 -0.955896854 -2.176721811 0.771060944 14.875996590 -28.476327896 -15.786120415 -5.331668218 22.538120270
0 -2.089187145 1.736116767 1.907819748 0 0 -1 -0.717318892 -0.734904408 0.904034972 -28.532371521 13.714955330 0.387423307 -7.276362737 11.813765526
"""
        
        with open('/tmp/test_rumd.xyz', 'w') as fh:
            fh.write(ioformat_1)
            self.input_file = fh.name

    def test_read_write(self):
        with TrajectoryRUMD(self.input_file) as th:
            with TrajectoryRUMD(self.input_file + 'out', 'w') as th_out:
                self.assertEqual(th.timestep, 0.005)
                self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
                th_out.timestep = th.timestep
                for i, system in enumerate(th):
                    th_out.write(system, th.steps[i])
        with TrajectoryRUMD(self.input_file + 'out') as th:
            self.assertEqual(th.timestep, 0.005)
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])


class TestUtils(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):
        with open('/tmp/test_1.xyz', 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:1 cell:5.0,5.0,5.0
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:2 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:4 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")
        with open('/tmp/test_2.xyz', 'w') as fh:
            fh.write("""\
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:2 cell:5.0,5.0,5.0
B 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:3 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:4 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
2
metafmt:space,comma columns:id,x,y,z mass:1.0,2.0,3.0 step:5 cell:5.0,5.0,5.0 
C 1.0 -1.0 0.0
B 2.9 -2.9 0.0
""")

    def test_paste(self):
        import atooms.trajectory as trj
        t1 = trj.TrajectoryXYZ('/tmp/test_1.xyz')
        t2 = trj.TrajectoryXYZ('/tmp/test_2.xyz')
        steps = []
        for step, s1, s2 in trj.utils.paste(t1, t2):
            steps.append(step)
        self.assertEqual(steps, [2, 4])

if __name__ == '__main__':
    unittest.main()


