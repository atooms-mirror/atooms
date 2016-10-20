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

    def test_no_output(self):
        # Disable writers completely
        s = Simulation(output_path=None, steps = 10, enable_speedometer=False)
        s.run()
        self.assertEqual(len(s._non_targeters),0)

    def test_target(self):
        s = Simulation(output_path='/tmp/test_simulation', steps = 10000, thermo_interval = 2000, config_number = 10)
        s.run()

    def test_target_multi(self):
        s = Simulation(output_path='/tmp/test_simulation', steps=10000, thermo_interval=2000, config_number=10)
        s.run()
        s.run(20000)

    def test_target_restart(self):
        s=Simulation(output_path='/tmp/test_simulation_restart', thermo_interval=2000, config_number=10)
        s.restart=True
        s.target_steps=10000
        s.run()
        s.target_steps=20000
        s.run()

    def test_target_restart_fake(self):
        s=Simulation(output_path='/tmp/test_simulation_restart', thermo_interval=2000, config_number=10)
        s.restart=True
        s.target_steps=10000
        s.run()
        s.target_steps=10000
        s.run()

if __name__ == '__main__':
    unittest.main()


