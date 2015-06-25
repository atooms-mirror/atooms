#!/usr/bin/env python

import unittest
from atooms.voronoi import TrajectoryVoronoi

class TestVoronoi(unittest.TestCase):

    def setUp(self):
        pass

    def test_read(self):
        t = TrajectoryVoronoi('reference/wahn_voronoi.xyz')
        v0 = []
        for i in t.samples:
            p, v = t.read_sample(i)
            v0.append(v[0].signature)
        ref = [[0, 1, 10, 4], [0, 3, 6, 7], [1, 1, 8, 3, 1], [0, 2, 8, 4], [1, 5, 2, 5, 3], [0, 1, 10, 4], [0, 2, 8, 4], [0, 3, 7, 5, 1], [2, 1, 7, 4, 1, 1], [0, 1, 10, 4], [0, 3, 6, 4], [0, 2, 8, 5], [0, 2, 8, 6], [2, 3, 3, 5, 1, 1], [0, 2, 8, 5], [0, 2, 8, 5], [0, 1, 10, 4], [4, 3, 2, 3, 2, 3], [0, 1, 10, 3], [1, 1, 8, 3, 1]]
        self.assertEqual(ref, v0)

if __name__ == '__main__':
    unittest.main(verbosity=0)


