#!/usr/bin/env python

import unittest
import numpy
from atooms import trajectory, postprocessing

def deviation(x, y):
    return (numpy.sum((x-y)**2)/len(x))**0.5

class TestMeanSquareDisplacement(unittest.TestCase):

    # TODO: speed up this test
    def test_msd(self):
        f = 'reference/kalj-matrix.h5'
        t = trajectory.Sliced(trajectory.TrajectoryHDF5(f), slice(0,1000,1))
        p = postprocessing.MeanSquareDisplacement(t, [0.0, 1.0, 10.0, 100.0])
        p.do()
        ref_grid = numpy.array([0, 0.992, 9.976, 100])
        ref_value = numpy.array([0.0, 0.080864, 0.325229, 2.5653])
        self.assertLess(deviation(p.grid, ref_grid), 1e-2)
        self.assertLess(deviation(p.value, ref_value), 1e-2)

class TestFourierSpace(unittest.TestCase):

    # TODO: this seed reset is not really working isnt it?
    def setUp(self):
        numpy.random.seed(10)

    # TODO: speed up this test
    def test_fskt(self):
        t = trajectory.TrajectoryHDF5('reference/lj.h5')
        p = postprocessing.SelfIntermediateScattering(t, kgrid=[7.0], tgrid=[0.0, 0.1, 0.5])
        p.do()
        ref_value = numpy.array([1.0, 0.51126058513090678, 0.017393617074980577])
        # Make sure p.value is numpy array
        self.assertLess(deviation(numpy.array(p.value[0]), ref_value), 0.01)

    # TODO: speed up this test
    def test_fskt_elongated(self):
        t = trajectory.TrajectoryHDF5('reference/lj_elongated/config.dat')
        p = postprocessing.SelfIntermediateScattering(t, kgrid=[7.0], tgrid=[0.0, 0.1, 0.5])
        p.do()
        ref_value = numpy.array([1.0, 0.70889408516023678, 0.18584564453072067])
        self.assertLess(deviation(numpy.array(p.value[0]), ref_value), 0.01)

class TestOverlap(unittest.TestCase):

    def test_overlap_pinned(self):
        t = trajectory.TrajectoryHDF5('reference/kalj_rumd_pinned.h5')
        p = postprocessing.OverlapDistribution(t, [20], skip=100)
        p.do()
        p.grid
        p.value


if __name__ == '__main__':
    unittest.main()


