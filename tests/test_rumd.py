#!/usr/bin/env python

import sys
import os
import unittest
from atooms.backends.rumd_backend import RumdBackend, single, multi
from atooms.simulation import Simulation, log
from atooms.simulation.parallel_tempering import ParallelTempering

log.setLevel(40)

def potential():
    import rumd
    pot = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
    pot.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
    pot.SetParams(i=1, j=0, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=0, j=1, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
    pot.SetParams(i=1, j=1, Epsilon=0.5, Sigma=0.88, Rcut=2.5)
    return [pot]

class Test(unittest.TestCase):

    def setUp(self):
        self.input_file = os.path.join(os.path.dirname(sys.argv[0]), '../data/KA_N256_L6.0_T0.80.xyz.gz')
            
    def test_single(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s), output_path='/tmp/test_rumd_single', steps=2000, 
                        thermo_interval=100, config_interval=100, checkpoint_interval=100,
                        restart=False)
        #si = Simulation(s, output_path=None, steps=200000) 
        si.run()

    def test_multi(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s), output_path='/tmp/test_rumd_multi', steps=2000, 
                        thermo_interval=100, config_interval=100, checkpoint_interval=100,
                        restart=False)
        si.run()
        si.run(1000)

    def test_multi_2(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s), output_path='/tmp/test_rumd_multi', steps=2000, 
                        thermo_interval=100, config_interval=100, checkpoint_interval=100,
                        restart=False)
        si.run()
        s = single(s, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s), output_path='/tmp/test_rumd_multi', steps=1000, 
                        thermo_interval=100, config_interval=100, checkpoint_interval=100,
                        restart=False)
        si.run()

    def test_multi_rmsd(self):
        s = single(self.input_file, potential, T=0.80, dt=0.002)
        si = Simulation(RumdBackend(s), output_path='/tmp/test_rumd_multi', 
                        thermo_interval=100, config_interval=100, checkpoint_interval=100, steps=1000000000,
                        restart=False)
        si.run(rmsd=1.0)
        si.run(rmsd=2.0)

    def test_pt(self):
        T = [1.0,0.95,0.9]
        nr = len(T)
        sa = [Simulation(RumdBackend(s, fixcm_interval=1000)) for s in multi(
            [self.input_file]*nr, potential, T=T, dt=[0.002]*nr)]
        pt = ParallelTempering(sa, T, '/tmp/test_rumd_rx', steps=5,
                               swap_interval=10000, thermo_interval=1, config_interval=1, 
                               checkpoint_interval=1, mute_config_except=[0], restart=False)
        pt.run()

    def tearDown(self):
        os.system('rm -rf /tmp/test_rumd_single')
        os.system('rm -rf /tmp/test_rumd_rx')

if __name__ == '__main__':
    unittest.main()
