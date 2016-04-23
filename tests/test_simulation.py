#!/usr/bin/env python

import unittest
from atooms.simulation import Simulation
from atooms.adapters.backend import SimulationBackend


class TestSimulationBackend(unittest.TestCase):

    def setUp(self):
        pass

    def test_simulation(self):
        s = Simulation(SimulationBackend(), None, steps = 10000)
        s.verbosity = 2
        s.run()

    def test_target(self):
        s = Simulation(SimulationBackend(), None, steps = 10000, thermo_interval = 20, config_number = 10)
        s.run()

if __name__ == '__main__':
    unittest.main()


