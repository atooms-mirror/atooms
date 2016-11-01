#!/usr/bin/env python

"""Test trajectory decorators."""

import unittest
import numpy

from atooms.trajectory.decorators import Unfolded
from atooms.trajectory import TrajectoryXYZ

class TestDecorators(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):        
        self.finp = '/tmp/test_decorators.xyz'
        with open(self.finp, 'w') as fh:
            fh.write("""\
2
step:1
A 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
step:2
A 1.1 -1.1 0.0
A -2.9 -2.9 0.0
2
step:3
A 1.2 -1.2 0.0
A -2.9 2.9 0.0
6.0 6.0 6.0
""")

    def test_unfolded_jump(self):
        with self.Trajectory(self.finp) as th:
            with Unfolded(th) as tu:
                self.assertEqual(list(tu[2].particle[1].position), [3.1, -3.1, 0.0]) # unfolded
                self.assertEqual(list(th[2].particle[1].position), [-2.9, 2.9, 0.0]) # original

        
if __name__ == '__main__':
    unittest.main(verbosity=0)


