#!/usr/bin/env python

import os
import sys
import unittest

from atooms.adapters.atoomsf90 import Simulation, System

# TODO: refactor adapters tests, they should be unique (except for constructors) because the interface is the same

# TODO: skip this if atoomsf90 is not there

class TestAtooms(unittest.TestCase):

    def setUp(self):
        self.file_ref = 'reference/kalj/config.dat'
        self.ref_u = -4805.69761109

    def test_potential_energy(self):
        s = System(self.file_ref)
        self.assertLess(abs(s.potential_energy()-self.ref_u), 1e-5)

    def test_potential_energy_from_simulation(self):
        s = Simulation(self.file_ref)
        self.assertLess(abs(s.system.potential_energy()-self.ref_u), 1e-5)

    def test_simulation_run(self):
        ref = -4740.46362485
        file_output = '/tmp/test_atooms_run.h5'
        s = Simulation(file_output, file_input=self.file_ref, 
                       opts={'--dt':0.001, '-c':1, '-e':1})
        s.setup(target_steps=10)
        s.run()
        self.assertLess(abs(s.system.potential_energy()-ref), 1e-5)
        os.remove(file_output)

if __name__ == '__main__':
    unittest.main(verbosity=0)
