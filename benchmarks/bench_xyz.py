#!/usr/bin/env python

import benchmark
import numpy

from atooms import trajectory 

class BenchmarkRead(benchmark.Benchmark):

    each = 3

    def setUp(self):
        pass

    def test_xyz_indexed(self):
        f = 'atooms/reference/gcm.xyz'
        t = trajectory.TrajectoryXYZIndexed(f)
        t.read_initial_state()
        for i in t.samples:
            t.read_sample(i)
        t.close()

    def test_xyz(self):
        f = 'atooms/reference/gcm.xyz'
        t = trajectory.TrajectoryXYZ(f)
        t.read_initial_state()
        for i in t.samples:
            t.read_sample(i)
        t.close()

class BenchmarkUnfold(benchmark.Benchmark):

    each = 3

    def setUp(self):
        pass

    def test_xyz_indexed(self):
        f = 'atooms/reference/gcm.xyz'
        t = trajectory.Unfolded(trajectory.TrajectoryXYZIndexed(f))
        t.read_initial_state()
        for i in t.samples:
            t.read_sample(i)
        t.close()

    def test_xyz(self):
        f = 'atooms/reference/gcm.xyz'
        t = trajectory.TrajectoryXYZ(f)
        t.read_initial_state()
        for i in t.samples:
            t.read_sample(i, unfolded=True)
        t.close()

if __name__ == '__main__':
    benchmark.main(format="markdown", numberFormat="%.4g")


