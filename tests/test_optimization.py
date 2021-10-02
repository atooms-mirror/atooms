#!/usr/bin/env python

import unittest

from atooms.optimization import Optimization
import atooms.backends.lammps
from atooms.backends.lammps import EnergyMinimization

SKIP = not atooms.backends.lammps.installed()

class Test(unittest.TestCase):

    def setUp(self):
        import os
        if SKIP:
            self.skipTest('missing LAMMPS')
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/lj_N256_rho1.0.xyz')

    def test_min(self):
        cmd = """
        pair_style      lj/cut 2.5
        pair_modify     shift yes
        pair_coeff      1 1 1.0 1.0 2.5
        """
        bck = EnergyMinimization(self.input_file, cmd)
        bck.verbose = False
        opt = Optimization(bck)
        opt.tolerance = 1e-4
        e = bck.system.potential_energy(per_particle=True)
        opt.run()
        e = bck.system.potential_energy(per_particle=True)
        w = bck.system.force_norm_square(per_particle=True)
        self.assertLess(e, -6.8)
        self.assertLess(w, 1e-4)


if __name__ == '__main__':
    unittest.main()
