#!/usr/bin/env python

import unittest
import numpy
from atooms.trajectory.field import TrajectoryField

class TestField(unittest.TestCase):

    def setUp(self):
        pass

    def test_write(self):
        with TrajectoryField('/tmp/test_field.xyz', 'w') as th:
            th.write([1.0, 2.0], step=0)
            th.write([-1.0, -2.0], step=1)
        with TrajectoryField('/tmp/test_field.xyz') as th:
            self.assertAlmostEqual(th[0][0], 1.)
#            self.assertEqual(th[1], numpy.array([-1., -2.]))

if __name__ == '__main__':
    unittest.main()


