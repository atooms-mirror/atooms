#!/usr/bin/env python

import unittest

from atooms.optimization import Optimization
try:
    from atooms.backends.lammps import EnergyMinimization
    SKIP = False
except ImportError:
    SKIP = True

class Test(unittest.TestCase):

    def setUp(self):
        import os
        if SKIP:
            self.skipTest('missing LAMMPS')
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/lj_N1000_rho1.0.xyz')

    def test_min(self):
        cmd = """
        pair_style      lj/cut 2.5
        pair_modify     shift yes
        pair_coeff      1 1 1.0 1.0 2.5
        """
        bck = EnergyMinimization(self.input_file, cmd)
        bck.verbose = False
        opt = Optimization(bck)
        e = bck.system.potential_energy(normed=True)
        #print('Initial energy {}'.format(e))
        opt.run()
        e = bck.system.potential_energy(normed=True)
        #print('Final energy {}'.format(e))
        #self.assertLessThan(e, 1e-10)

if __name__ == '__main__':
    unittest.main()


