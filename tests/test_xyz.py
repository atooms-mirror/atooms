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
        t = trajectory.TrajectoryXYZIndexed(self.finp)
        self.assertEqual(t.samples, [0, 1, 2, 3])
        self.assertEqual(t.steps, [1, 2, 3, 4])
        self.assertEqual(r_ref[0], list(t.read_sample(0).particle[0].position))
        self.assertEqual(r_ref[1], list(t.read_sample(0).particle[1].position))

    def test_xyz_indexed_unfolded(self):
        t1 = trajectory.TrajectoryXYZIndexed(self.finp)
        t = trajectory.Unfolded(t1)
        t.read_initial_state()
        t.read_sample(0)
        t.read_sample(1)
        self.assertEqual(list(t.read_sample(2).particle[1].position), [3.1, -3.1, 0.0])
        self.assertEqual(list(t1.read_sample(2).particle[1].position), [-2.9, 2.9, 0.0])
        t1.close()
        t.close()

    def test_xyz_indexed_unfolded_skip(self):
        """Test that unfolded trajectories can skip samples"""
        t1 = trajectory.TrajectoryXYZIndexed(self.finp)

        # Here we read all samples
        t = trajectory.Unfolded(t1)
        t.read_initial_state()
        t.read_sample(0)
        t.read_sample(1)
        t.read_sample(2)
        s = t.read_sample(3)

        # Here we skip one sample
        t1 = trajectory.Unfolded(t1)
        t1.read_initial_state()
        t1.read_sample(0)
        s1 = t1.read_sample(3)

        self.assertEqual(list(s.particle[0].position), list(s1.particle[0].position))
        self.assertEqual(list(s.particle[1].position), list(s1.particle[1].position))

        t1.close()
        t.close()

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
        #self.assertEqual(list(t.read_sample(1).particle[0].position), [1.1, 1.2, 1.3])
        #self.assertEqual(list(t.read_sample(1).particle[1].position), [-1.1, -1.2, -1.3])
        
if __name__ == '__main__':
    unittest.main(verbosity=0)


