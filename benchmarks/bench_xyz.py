#!/usr/bin/env python

import benchmark
import numpy

from atooms import trajectory 

class BenchmarkRead(benchmark.Benchmark):

    each = 1

    def setUp(self):
        pass

    def test_xyz(self):
        f = 'reference/gcm.xyz'
        f = 'data/trajectory.xyz'
        with trajectory.TrajectoryXYZ(f) as t:
            for s in t[:40]:
                pass

    def test_simple_xyz(self):
        f = 'reference/gcm.xyz'
        f = 'data/trajectory.xyz'
        with trajectory.TrajectorySimpleXYZ(f) as t:
            for s in t[:40]:
                pass

# class BenchmarkUnfold(benchmark.Benchmark):

#     each = 3

#     def setUp(self):
#         pass

#     def test_xyz_indexed(self):
#         f = 'reference/gcm.xyz'
#         with trajectory.Unfolded(trajectory.TrajectoryXYZ(f)) as t:
#             for s in t:
#                 pass


if __name__ == '__main__':
    benchmark.main(format="markdown", numberFormat="%.4g")


