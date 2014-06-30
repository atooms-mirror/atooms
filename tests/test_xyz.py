#!/usr/bin/env python

import unittest
import numpy

from atooms import trajectory 

class TestXYZ(unittest.TestCase):

    def setUp(self):
        self.finp = '/tmp/test_pbc.xyz'
        fh = open(self.finp, 'w')
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
        fh.close()

    def test_xyz_indexed(self):
        r_ref = [[1., -1., 0.], [2.9, -2.9, 0.]]
        t = trajectory.TrajectoryXYZ(self.finp)
        self.assertEqual(t.steps, [1, 2, 3, 4])
        self.assertEqual(r_ref[0], list(t[0].particle[0].position))
        self.assertEqual(r_ref[1], list(t[0].particle[1].position))

    def test_xyz_indexed_unfolded(self):
        t1 = trajectory.TrajectoryXYZ(self.finp)
        t = trajectory.Unfolded(t1)
        t[0]
        t[1]
        self.assertEqual(list(t[2].particle[1].position), [3.1, -3.1, 0.0])
        self.assertEqual(list(t1[2].particle[1].position), [-2.9, 2.9, 0.0])
        t1.close()
        t.close()

    def test_xyz_indexed_unfolded_skip(self):
        """Test that unfolded trajectories can skip samples"""
        # Here we read all samples
        t1 = trajectory.TrajectoryXYZ(self.finp)
        t = trajectory.Unfolded(t1)
        t[0]
        t[1]
        t[2]
        s = t[3]
        t1.close()
        t.close()

        # Here we skip one sample
        t1 = trajectory.TrajectoryXYZ(self.finp)
        t1 = trajectory.Unfolded(t1)
        t1[0]
        s1 = t1[3]
        t1.close()
        t.close()

        self.assertEqual(list(s.particle[0].position), list(s1.particle[0].position))
        self.assertEqual(list(s.particle[1].position), list(s1.particle[1].position))


    @unittest.skip('ikeda2 does not work')
    def test_ikeda2_indexed(self):
        f = '/tmp/test_xyz_ikeda2.xyz'
        with open('/tmp/test_xyz_ikeda2.xyz', 'w') as fh:
            fh.write("""\
           2            2    0.12000000E+02    0.50000000E+03    0.20000000E+01
    0.50050000E+07
               1    0.1 0.2 0.3 1.1 1.2 1.3
               2    -0.1 -0.2 -0.3 1.1 1.2 1.3
    0.50100000E+07
               1    1.1 1.2 1.3 1.1 1.2 1.3
               2    -1.1 -1.2 -1.3 1.1 1.2 1.3
""")
        t = trajectory.TrajectoryXYZIkeda2Indexed(f)
        self.assertEqual(t.steps, [5005000, 5010000])

    def test_xyz_with(self):
        r_ref = [[1., -1., 0.], [2.9, -2.9, 0.]]
        with trajectory.TrajectoryXYZ(self.finp) as t:
            self.assertEqual(r_ref[0], list(t[0].particle[0].position))
            self.assertEqual(r_ref[1], list(t[0].particle[1].position))

        with trajectory.TrajectoryXYZ(self.finp + '.out', 'w') as to:
            with trajectory.TrajectoryXYZ(self.finp) as t:
                to.write(t[0], 0)
                
        with trajectory.TrajectoryXYZ(self.finp + '.out') as to:
            self.assertEqual(r_ref[0], list(to[0].particle[0].position))

    def test_xyz_iter(self):
        with trajectory.TrajectoryXYZ(self.finp) as t:
            for s in t:
                s.particle
            t[-1]
        
if __name__ == '__main__':
    unittest.main(verbosity=0)


