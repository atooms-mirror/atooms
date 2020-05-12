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
from atooms.trajectory.rumd import TrajectoryRUMD, SuperTrajectoryRUMD

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
        from atooms.backends.rumd import System
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        s = System(self.backend.rumd_simulation.sample)
        si.add(write_thermo, 100)
        si.add(write_config, 100)
        si.run()
        ls = glob.glob('/tmp/test_rumd_single/trajectory/*')
        self.assertEqual(len(ls), 21)

    @unittest.expectedFailure
    def test_species(self):
        """
        This test is known to fail because RUMD sample do not return the
        particle types as array and therefore they are not propagated
        when particles are not sorted by species
        """
        from atooms.backends.rumd import System
        # Create a new input file with one particle species changed
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/ka_N256_rho1.185_rumd.xyz.gz')
        from atooms.core.utils import mkdir
        mkdir('/tmp/test_rumd_species/')
        with TrajectoryRUMD(self.input_file) as th:
            system = th[-1]
        system.particle[0].species = system.particle[-1].species
        with TrajectoryRUMD('/tmp/test_rumd_species/input.xyz.gz', 'w') as th:
            th.write(system)

        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_species/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        self.assertEqual(system.particle[0].species, system.particle[-1].species)
        si.add(write_config, 100)
        si.run()
        with SuperTrajectoryRUMD('/tmp/test_rumd_species/trajectory') as th:
            system = th[-1]
        self.assertEqual(system.particle[0].species, system.particle[-1].species, 'rumd requires ordered types')

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

    def test_system_copy(self):
        """Make sure rumd systems can be deepcopied"""
        import copy
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        s0 = copy.copy(si.system)
        s1 = copy.deepcopy(si.system)
        self.assertEqual(si.system.particle[0].position[-1], s1.particle[0].position[-1])

    def test_potential(self):
        import copy
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        pos0 = si.system.particle[0].position[0]
        s = copy.deepcopy(si.system)
        si.run(100)
        pos1 = si.system.particle[0].position[0]
        si.system = s
        si.run(100)
        pos2 = si.system.particle[0].position[0]
        self.assertTrue(abs(pos1 - pos0)>1e-2)
        self.assertTrue(abs(pos1 - pos2)<1e-4)

    def test_copy_and_run(self):
        self.skipTest('skipped test')
        import copy
        import atooms.trajectory.ram
        from atooms.backends.rumd import RUMD
        import rumd
        potential = rumd.Pot_LJ_12_6(cutoff_method=rumd.ShiftedPotential)
        potential.SetVerbose(False)
        potential.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
        potential.SetParams(i=1, j=0, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
        potential.SetParams(i=0, j=1, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
        potential.SetParams(i=1, j=1, Epsilon=0.5, Sigma=0.88, Rcut=2.5)
        input_file = os.path.join(os.path.dirname(__file__),
                                  '../data/ka_N256_rho1.185_rumd.xyz.gz')
        backend = RUMD(input_file, integrator='nvt', temperature=0.8, dt=0.002, forcefield=[potential])
        sim = Simulation(backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=100, restart=False)

        # First set if runs
        # print sim.system.particle[0].position[0]
        tr = atooms.trajectory.ram.TrajectoryRamFull()
        th = atooms.trajectory.ram.TrajectoryRamFull()
        tr[0] = sim.system
        # print tr[0].particle[0].position[0]
        sim.run(1000)  # sometimes is nan, depending on steps

        # Ram does not work with rumd because of change in System init
        th[0] = sim.system
        # print sim.system.particle[0].position, sim.system.potential_energy(per_particle=True)

        sim.run(100)
        # print sim.system.particle[0].position[0]
        # this assignment leads to trouble, BUT ONLY IF WE STORE THE SAMPLES IN tr TRAJECTORY
        # assigning tr[0] above also leads to a different trajectory...!!
        sim.system = th[0]
        #sim.system = tr[0]
        # print sim.system.particle[0].position, sim.system.potential_energy(per_particle=True)

    def test_unfold(self):
        from atooms.backends.rumd import unfold
        def unf(sim):
            unfold(sim.system).particle[0].position
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, restart=False)
        si.add(unf, 100)
        si.run()

    def test_leakage(self):
        #self.skipTest('skipped test')
        from atooms.backends.rumd import System
        from atooms.trajectory.ram import TrajectoryRamFull
        from atooms.backends.rumd import unfold
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        # self.backend.rumd_simulation.sample.__swig_destroy__(self.backend.rumd_simulation.sample)
        #del self.backend.rumd_simulation.sample

        # This does not leak memory
        trj = TrajectoryRamFull()
        trj[0] = si.system
        for i in range(5):
            si.system = trj[0]
            si.run()
            unfold(si.system).particle[0].position[0], unfold(trj[0]).particle[0].position[0]

        # This does not leak memory anymore because we use System.update()
        # trj = TrajectoryRamFull()
        # for i in range(5):
        #     trj[0] = si.system
        #     si.run()
        #     #print si.system.particle[0].position[0], trj[0].particle[0].position[0]
        #     print unfold(si.system).particle[0].position[0], unfold(trj[0]).particle[0].position[0]

    def test_update(self):
        #self.skipTest('skipped test')
        from atooms.backends.rumd import System
        from atooms.trajectory.ram import TrajectoryRamFull
        si = Simulation(self.backend,
                        output_path='/tmp/test_rumd_single/trajectory',
                        steps=2000, checkpoint_interval=100,
                        restart=False)
        itg = si.system.sample.GetIntegrator()
        trj = TrajectoryRamFull()
        trj[0] = si.system
        si.system.thermostat.temperature = 10.0
        # Objects are different but underlying information is shared...
        # print trj[0].sample.GetIntegrator().GetInfoString(18), itg.GetInfoString(18)
        self.assertNotEqual(trj[0].sample.GetIntegrator(), itg)
        #self.assertEqual(trj[0].sample.GetIntegrator().GetInfoString(18), itg.GetInfoString(18))
        trj[0] = si.system  # This uses update
        # print trj[0].sample.GetIntegrator().GetInfoString(18), itg.GetInfoString(18)
        self.assertNotEqual(trj[0].sample.GetIntegrator(), itg)

    def test_new_format(self):
        input_file = os.path.join(os.path.dirname(__file__),
                                  '../data/ka_rumd_N1200.xyz.gz')
        with TrajectoryRUMD(input_file) as th:
            self.assertAlmostEqual(th[0].cell.side[0], 10.03)
            self.assertAlmostEqual(th[0].cell.side[1], 10.03)
            self.assertAlmostEqual(th[0].cell.side[2], 10.03)

    def tearDown(self):
        os.system('rm -rf /tmp/test_rumd_*')


if __name__ == '__main__':
    unittest.main()
