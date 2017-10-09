#!/usr/bin/env python

import os
import unittest
from atooms.simulation import Simulation, write_thermo, write_config, target
from atooms.core.utils import setup_logging
try:
    from atooms.backends.lammps import LAMMPS, Interaction
    SKIP = False
except ImportError:
    SKIP = True

setup_logging(level=40)


class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing LAMMPS')
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/lj_N1000_rho1.0.xyz')

    def test_single(self):
        import sys
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        bck = LAMMPS(self.input_file, cmd)
        sim = Simulation(bck)
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.62635, places=5)
        sim.run(10)
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.64526, places=5)
        sim.run(10)
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.67598, places=5)

    def test_energy(self):
        import sys
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        """
        bck = LAMMPS(self.input_file, cmd)
        bck.system.interaction.compute("energy", bck.system.particle, bck.system.cell)
        self.assertEqual(bck.system.interaction.energy, -4.2446836)  # crosschecked

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
