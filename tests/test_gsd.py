#!/usr/bin/env python

import os
import unittest
import numpy as np

from atooms import trajectory


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test(self):
        fname = '../data/trajectory3d.gsd'
        traj = trajectory.TrajectoryGSD(fname)
        p = traj[0].particle[0]

        # Sanity checks.
        self.assertEqual(len(traj), 4)
        self.assertEqual(traj.steps, [0, 1, 2, 3])
        self.assertEqual(traj[0].distinct_species(), ['A', 'B'])
        self.assertTrue( np.all(np.isclose(traj[0].cell.side, np.array([10.68387318, 10.68387318, 10.68387318]), 1e-12)) )
        self.assertEqual(traj[0].number_of_dimensions, 3)
        
        self.assertEqual(p.species, 'B')
        self.assertEqual(p.mass, 1.0)
        self.assertEqual(p.radius, 0.5)
        self.assertTrue( np.all(p.position == np.array([5.3419366, 5.3419366, 5.3419366], dtype=np.float32)) )

    def test2d(self):
        fname = '../data/trajectory2d.gsd'
        traj = trajectory.TrajectoryGSD(fname)
        p = traj[0].particle[0]

        # Sanity checks.
        self.assertEqual(len(traj), 4)
        self.assertEqual(traj.steps, [0, 1, 2, 3])
        self.assertEqual(traj[0].distinct_species(), ['A', 'B'])
        self.assertTrue( np.all(np.isclose(traj[0].cell.side, np.array([10.78327751, 10.78327751]), 1e-12)) )
        self.assertEqual(traj[0].number_of_dimensions, 2)
        
        self.assertEqual(p.species, 'A')
        self.assertEqual(p.mass, 1.0)
        self.assertEqual(p.radius, 0.5)
        self.assertTrue( np.all(p.position == np.array([5.3916388, 5.3916388], dtype=np.float32)) )

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
