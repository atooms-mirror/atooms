#!/usr/bin/env python

import os
import glob
import unittest
try:
    import rumd
    SKIP = False
except ImportError:
    SKIP = True
from atooms.simulation import Simulation, Scheduler, WriterThermo, WriterConfig, TargetRMSD
from atooms.utils import setup_logging

setup_logging(level=40)

class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing RUMD')
        self.input_file = os.path.join(os.path.dirname(__file__), '../data/KA_N256_rumd.xyz.gz')
        self.forcefield_file = os.path.join(os.path.dirname(__file__), '../data/ka_rumd.ff')            
        from atooms.simulation.backends import RumdBackend
        self.backend = RumdBackend(self.input_file,
                                   self.forcefield_file, integrator='nvt',
                                   temperature=0.80, dt=0.002)
        self.backend_2 = RumdBackend(self.input_file,
                                     self.forcefield_file, integrator='nvt',
                                     temperature=0.80, dt=0.002)

    def test_single(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        ls = glob.glob('/tmp/test_rumd_single/trajectory/*')
        self.assertEqual(len(ls), 21)

    def test_multi(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        si.run(1000)

    def test_multi_writing(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_writing/trajectory',
                        enable_speedometer=False)
        si.add(WriterConfig(), Scheduler(500))
        si.run(50000)
        si.run(5000)
        ls = glob.glob('/tmp/test_rumd_multi_writing/trajectory/*')
        self.assertEqual(len(ls), 11)

    def test_multi_2(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        si = Simulation(self.backend_2,
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=1000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()

    def test_multi_rmsd(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_rmsd/trajectory',
                        checkpoint_interval=100, steps=1000000000,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.add(TargetRMSD(1.0), Scheduler(1000))
        si.run()

    def tearDown(self):
        os.system('rm -rf /tmp/test_rumd_*')

if __name__ == '__main__':
    unittest.main()
