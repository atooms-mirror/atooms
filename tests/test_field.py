#!/usr/bin/env python

import unittest
import numpy
from atooms.trajectory.xyz import TrajectoryXYZ
from atooms.trajectory.field import TrajectoryField

class TestField(unittest.TestCase):

    def setUp(self):
        with TrajectoryField('/tmp/test_field.xyz', 'w') as th:
            th.write([1.0, 2.0], step=0)
            th.write([-1.0, -2.0], step=1)

    def test_read(self):
        with TrajectoryField('/tmp/test_field.xyz') as th:
            self.assertAlmostEqual(th[0][0], 1.)
            self.assertAlmostEqual(th[1][1], -2.)

    def test_read_as_xyz(self):
        with TrajectoryXYZ('/tmp/test_field.xyz') as th:
            self.assertAlmostEqual(th[0].particle[0].field, 1.)
            self.assertAlmostEqual(th[1].particle[1].field, -2.)

if __name__ == '__main__':
    unittest.main()


