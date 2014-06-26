#!/usr/bin/env python

import unittest
from atooms import simulation


class TestSimulationBackend(unittest.TestCase):

    def setUp(self):
        from atooms.adapters import backend
        self.s = backend.Simulation()
        self.s.verbosity = 2

    def test_simulation(self):
        #self.s.add(simulation.TargetSteps(10000))
        self.s.setup(target_steps = 10000)
        self.s.run()

    def test_target(self):
        self.s.setup(target_steps = 100, thermo_period = 20, config_number = 10)
        self.s.run()

if __name__ == '__main__':
    unittest.main()


