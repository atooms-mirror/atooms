#!/usr/bin/env python

import fuzzyunittest
import numpy
from atooms.system.particle import Particle
from atooms.trajectory import TrajectoryNeighbors, Trajectory, TrajectoryHDF5
from atooms.postprocessing.boo import BondOrientationalOrder, periodic_vector


class TestBOO(fuzzyunittest.FuzzyTestCase):

    def setUp(self):
        pass

    def _load(self, structure):
        # TODO: these references files don't exist anymore
        if structure == 'bcc':
            self.t = Trajectory('reference/lj_bcc.h5')
            self.tn = TrajectoryNeighbors('reference/lj_bcc.h5.voronoi.xyz.neigh')
            self.ref = {6: 0.510688230857, 4: 0.0363696483727, 8: 0.429322472922}
        elif structure == 'fcc':
            self.t = Trajectory('reference/lj_fcc.h5.min')
            self.tn = TrajectoryNeighbors('reference/lj_fcc.h5.min.voronoi.xyz.neigh')
            self.ref = {6: 0.574524259714, 4: 0.190940653956, 8: 0.403914561085}
        elif structure == 'fluid':
            self.t = TrajectoryHDF5('reference/lj.h5')
            self.tn = TrajectoryNeighbors('reference/lj.h5.voronoi.xyz.neigh')
            self.ref = {6: 0.2731}

    def _test(self):
        s = self.t.read_initial_state()
        box = s.cell.side
        np = len(s.particle)
        for i in self.t.samples[:1]:
            s = self.t.read_sample(i)
            n = self.tn.read_sample(i)
            boo = BondOrientationalOrder(s.particle, n, box)
            #r = numpy.array([p.position for p in s.particle])
            #boo = BondOrientationalOrder(r, numpy.array(n), box)
            # import copy
            # for pi, q, j in zip(s.particle, boo.ql(6), n): # / np, numpy.sum(boo.ql(4)) / np
            #     pn = copy.deepcopy([s.particle[jj] for jj in j])
            #     ri = pi.position
            #     rkall = [pk.nearest_image(pi, s.cell).position for pk in pn]
            #     dr = [numpy.sqrt(numpy.dot(ri-rk, ri-rk)) for rk in rkall]
                #print q , dr
            print numpy.average(boo.ql(6))
        for l in self.ref:
            self.assertAlmostEqual(self.ref[l], numpy.average(boo.ql(l)))

    # def test_fcc(self):
    #     self._load('fcc')
    #     self._test()

    # def test_bcc(self):
    #     self._load('bcc')
    #     self._test()

    def test_fluid(self):
        self._load('fluid')
        self._test()
        

if __name__ == '__main__':
    unittest.main(verbosity=0)


