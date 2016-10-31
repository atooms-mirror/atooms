#!/usr/bin/env python

import unittest
import numpy

from atooms.trajectory import Unfolded
from atooms.trajectory import TrajectoryXYZ, TrajectorySimpleXYZ, TrajectoryNewXYZ

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

class TestSimpleXYZ(TestXYZ):

    # TODO: refactor generic tests for trajectories like this :-)
    Trajectory = TrajectorySimpleXYZ

    def test_xyz_meta(self):
        pass

    def test_xyz_mass(self):
        pass

        
if __name__ == '__main__':
    unittest.main(verbosity=0)


