#!/usr/bin/env python

import unittest
import logging
import numpy
from atooms.simulation import Simulation, WriterThermo, Scheduler, write_thermo
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
        s = Simulation(output_path='/tmp/test_simulation/trajectory', steps=100)
        s.run()

    def test_target_restart(self):
        f='/tmp/test_simulation_restart/trajectory'
        s=Simulation(output_path=f)
        s.add(write_thermo, Scheduler(20))
        s.run(100)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 100)

        s=Simulation(output_path=f, restart=True)
        s.add(write_thermo, Scheduler(20))
        s.run(200)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 200)

    def test_target_restart_fake(self):
        f = '/tmp/test_simulation_restart/trajectory'
        s=Simulation(output_path=f)
        #s.add(WriterThermo(), Scheduler(20))
        s.add(write_thermo, Scheduler(20))
        s.run(100)
        s.run(100)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 100)

    def test_scheduler(self):
        class Simulation:
            def __init__(self):
                self.steps = 0
        s = Scheduler(3)
        sim = Simulation()
        inext = []
        for i in range(8):
            sim.steps = i
            inext.append(s(sim))

        self.assertEqual(inext, [3, 3, 3, 6, 6, 6, 9, 9])

if __name__ == '__main__':
    unittest.main()


