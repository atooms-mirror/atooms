#!/usr/bin/env python

import unittest
import logging
import numpy
from atooms.simulation import Simulation
from atooms.utils import setup_logging

setup_logging(level=40)

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_no_output(self):
        # Disable writers completely
        s = Simulation(output_path=None, steps=10, enable_speedometer=False)
        s.run()
        self.assertEqual(len(s._non_targeters),0)

    def test_target(self):
        s = Simulation(output_path='/tmp/test_simulation/trajectory', steps=100, thermo_interval=20)
        s.run()

    def test_target_restart(self):
        f='/tmp/test_simulation_restart/trajectory'
        s=Simulation(output_path=f, thermo_interval=20)
        s.target_steps=100
        s.run()
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 100)
        s=Simulation(output_path=f, thermo_interval=20, restart=True)
        s.run(200)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 200)

    def test_target_restart_fake(self):
        f = '/tmp/test_simulation_restart/trajectory'
        s=Simulation(output_path=f, thermo_interval=20,)
        s.target_steps=100
        s.run()
        s.target_steps=100
        s.run()
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 100)

if __name__ == '__main__':
    unittest.main()


