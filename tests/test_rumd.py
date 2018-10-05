#!/usr/bin/env python

import os
import glob
import unittest
try:
    import rumd
    SKIP = False
except ImportError:
    SKIP = True
from atooms.simulation import Simulation, write_thermo, write_config, target
from atooms.core.utils import setup_logging

setup_logging(level=40)


class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing RUMD')
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/ka_N256_rho1.185_rumd.xyz.gz')
        self.forcefield_file = os.path.join(os.path.dirname(__file__),
                                            '../data/ka_rumd.ff')
        from atooms.backends.rumd import RUMD
        self.backend = RUMD(self.input_file,
                            self.forcefield_file, integrator='nvt',
                            temperature=0.80, dt=0.002)
        self.backend_2 = RUMD(self.input_file,
                              self.forcefield_file, integrator='nvt',
                              temperature=0.80, dt=0.002)

    def test_properties(self):
        from atooms.trajectory import TrajectoryRUMD
        t = TrajectoryRUMD(self.input_file)
        s0 = t[-1]
        sim = Simulation(self.backend,
                       output_path='/tmp/test_rumd_temp/trajectory',
                       steps=1, restart=False)
        s1 = sim.system
        self.assertAlmostEqual(s1.temperature, s0.temperature)
        self.assertAlmostEqual(s1.cell.side[0], s0.cell.side[0])
        self.assertAlmostEqual(s1.particle[0].position[0], s0.particle[0].position[0])

    def test_single(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(write_thermo, 100)
        si.add(write_config, 100)
        si.run()
        ls = glob.glob('/tmp/test_rumd_single/trajectory/*')
        self.assertEqual(len(ls), 21)

    def test_multi(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(write_thermo, 50)
        si.add(write_config, 100)
        si.run()
        si.run(1000)
        ls = glob.glob('/tmp/test_rumd_multi/trajectory/*')
        self.assertEqual(len(ls), 31)
        self.assertEqual(si.current_step, 3000)
        tmp = open('/tmp/test_rumd_multi/trajectory.thermo', 'r').readlines()
        steps = int(tmp[-1].split()[0])
        self.assertEqual(steps, 3000)
        self.assertEqual(len(tmp), 61 + 1)  # one is for comment line

    def test_multi_writing(self):
        # Test that we cumulate current_step and configurations
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_writing/trajectory',
                        enable_speedometer=False)
        si.add(write_config, 500)
        si.run(10000)
        si.run(5000)
        ls = glob.glob('/tmp/test_rumd_multi_writing/trajectory/*')
        self.assertEqual(len(ls), 31)
        self.assertEqual(si.current_step, 15000)

        # This second simulation tests that the folder gets cleared up
        si = Simulation(self.backend, steps=5000,
                        output_path='/tmp/test_rumd_multi_writing/trajectory',
                        enable_speedometer=False)
        si.add(write_config, 500)
        si.run()
        si.run()
        ls = glob.glob('/tmp/test_rumd_multi_writing/trajectory/*')
        self.assertEqual(len(ls), 21)
        self.assertEqual(si.current_step, 10000)

    def test_multi_2(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        si.add(write_thermo, 100)
        si.add(write_config, 100)
        si.run()
        si = Simulation(self.backend_2,
                        output_path='/tmp/test_rumd_multi_2/trajectory',
                        steps=1000, checkpoint_interval=100,
                        restart=False)
        si.add(write_thermo, 100)
        si.add(write_config, 100)
        si.run()

    def test_multi_rmsd(self):
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_multi_rmsd/trajectory',
                        checkpoint_interval=100, steps=1000000000,
                        restart=False)
        si.add(write_thermo, 100)
        si.add(write_config, 100)
        si.add(target, 1000, 'rmsd', 1.0)
        si.run()

    def tearDown(self):
        os.system('rm -rf /tmp/test_rumd_*')


if __name__ == '__main__':
    unittest.main()
