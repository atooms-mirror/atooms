#!/usr/bin/env python

import unittest
from atooms import simulation


class TestSimulationBackend(unittest.TestCase):

    def setUp(self):
        pass

    def test_simulation(self):
        from atooms.adapters import backend
        s = backend.Simulation(None, None, steps = 10000)
        s.verbosity = 2
        s.run()

    def test_target(self):
        from atooms.adapters import backend
        s = backend.Simulation(None, None, steps = 100, thermo_interval = 20, config_number = 10)
        s.run()

if __name__ == '__main__':
    unittest.main()


