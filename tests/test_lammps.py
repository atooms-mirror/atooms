#!/usr/bin/env python

import os
import unittest
from atooms.simulation import Simulation, write_thermo, write_config, target
from atooms.utils import setup_logging
try:
    from atooms.simulation.backends import LammpsBackend
    SKIP = False
except ImportError:
    SKIP = True

setup_logging(level=40)

class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing LAMMPS')
        self.input_file = os.path.join(os.path.dirname(__file__), '../data/lj_rho1.0.xyz')

    def test_single(self):
        import sys
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        bck = LammpsBackend(self.input_file, cmd)
        sim = Simulation(bck)
        sim.run(10)
        self.assertAlmostEqual(sim.system.particle[0].position[0], 3.64526, places=5)
        sim.run(10)
        self.assertAlmostEqual(sim.system.particle[0].position[0], 3.67598, places=5)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
