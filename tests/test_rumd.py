#!/usr/bin/env python

import os
import glob
import unittest
try:
    from atooms.backends.rumd_backend import RumdBackend, single, multi
    SKIP = False
except ImportError:
    SKIP = True
from atooms.simulation import Simulation, Scheduler, WriterThermo, WriterConfig, TargetRMSD
from atooms.utils import setup_logging

setup_logging(level=40)

def potential():
    import rumd
    pot = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
    pot.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
    pot.SetParams(i=1, j=0, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=0, j=1, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=1, j=1, Epsilon=0.5, Sigma=0.88, Rcut=2.5)
    pot.SetVerbose(False)
    return [pot]

class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing RUMD')
        self.input_file = os.path.join(os.path.dirname(__file__), '../data/KA_N256_rumd.xyz.gz')
            
    def test_single(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s),
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        ls = glob.glob('/tmp/test_rumd_single/trajectory/*')
        self.assertEqual(len(ls), 21)

    def test_multi(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s),
                        output_path='/tmp/test_rumd_multi/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        si.run(1000)

    def test_multi_writing(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002, interval_energy=500, interval_config=500)
        si = Simulation(RumdBackend(s),
                        output_path='/tmp/test_rumd_multi_writing/trajectory',
                        enable_speedometer=False)
        si.add(WriterConfig(), Scheduler(500))
        si.run(50000)
        si.run(5000)
        ls = glob.glob('/tmp/test_rumd_multi_writing/trajectory/*')
        self.assertEqual(len(ls), 11)

    def test_multi_2(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s),
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()
        s = single(s, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s),
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=1000, checkpoint_interval=100,
                        restart=False)
        si.add(WriterThermo(), Scheduler(100))
        si.add(WriterConfig(), Scheduler(100))
        si.run()

    def test_multi_rmsd(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s),
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
