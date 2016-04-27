#!/usr/bin/env python

import unittest
from atooms.simulation import Simulation, log

log.setLevel(40)

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_simulation(self):
        s = Simulation(steps = 10000)
        s.run()

    def test_target(self):
        s = Simulation(output_path='/tmp/test_simulation', steps = 10000, thermo_interval = 2000, config_number = 10)
        s.run()

    def test_target_multi(self):
        s = Simulation(output_path='/tmp/test_simulation', steps = 10000, thermo_interval = 2000, config_number = 10)
        s.run()
        s.run(20000)

if __name__ == '__main__':
    unittest.main()


