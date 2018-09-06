#!/usr/bin/env python

import unittest
import logging
import numpy
from atooms.simulation import Simulation, Scheduler, write_thermo, write_config, target_rmsd, write
from atooms.backends.dryrun import DryRun
from atooms.core.utils import setup_logging, rmd


setup_logging(level=40)

class Test(unittest.TestCase):

    def setUp(self):
        setup_logging(level=40, update=True)

    def test_no_output(self):
        # Disable writers completely
        s = Simulation(DryRun(), output_path=None, steps=10, enable_speedometer=False)
        s.run()
        self.assertEqual(len(s._non_targeters), 0)

    def test_multiple_run_calls(self):
        """
        Multiple calls to run() with varying number of steps should add up
        correctly. This was not the case in version <= 1.4.3.
        """
        from atooms.system import System

        # Minimal backend
        class Backend(object):

            def __init__(self):
                self.system = System()

            def run(self, steps):
                for i in range(steps):
                    pass

        backend = Backend()

        # The run_until() method should work correctly
        from atooms.simulation import Simulation
        simulation = Simulation(backend)
        simulation.run(10)
        simulation.run_until(30)
        self.assertEqual(simulation.current_step, 30)

        # The run() method called twice should also work correctly
        from atooms.simulation import Simulation
        simulation = Simulation(backend)
        simulation.run(10)
        simulation.run(20)
        self.assertEqual(simulation.current_step, 30)

    def test_target(self):
        s = Simulation(DryRun(), output_path='/tmp/test_simulation/trajectory', steps=100)
        s.run()

    def test_target_restart(self):
        f = '/tmp/test_simulation/restart/trajectory'
        s = Simulation(DryRun(), output_path=f)
        s.add(write_thermo, Scheduler(20))
        s.run(100)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 100)

        s = Simulation(DryRun(), output_path=f, restart=True)
        s.add(write_thermo, Scheduler(20))
        s.run(200)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 200)

    def test_config(self):
        from atooms.trajectory import TrajectoryXYZ
        f = '/tmp/test_simulation/config/trajectory.xyz'

        # We do not accept too deep introspection
        with self.assertRaises(ValueError):
            # Mute errors temporarily
            setup_logging(level=50, update=True)
            s = Simulation(DryRun(), output_path=f, enable_speedometer=False, steps=100)
            s.add(write, Scheduler(20), 'output', ['system.particle.position'])
            s.run()

        # Test generic writer and write_config
        setup_logging(level=40, update=True)
        s = Simulation(DryRun(), output_path=f, enable_speedometer=True, steps=100)
        s.trajectory = TrajectoryXYZ
        s.add(write_config, Scheduler(20))
        s.add(write, Scheduler(20), 'output', ['current_step',
                                               'system.cell'])
        s.run()
        import os
        self.assertTrue(os.path.exists(f))
        self.assertTrue(os.path.exists(f + '.output'))

    def test_target_rmsd(self):
        f = '/tmp/test_simulation/config/trajectory'
        with self.assertRaises(IndexError):
            s = Simulation(DryRun(), output_path=f, steps=100)
            s.add(target_rmsd, Scheduler(20))
            s.run()
        s = Simulation(DryRun(), output_path=f, steps=100)
        s.add(target_rmsd, Scheduler(20), 1.0)
        s.run()

    def test_target_walltime(self):
        """Check that walltime targeting works."""
        from atooms.simulation.observers import target_walltime
        f = '/tmp/test_simulation/config/trajectory'
        s = Simulation(DryRun(), output_path=f, steps=1000000000)
        s.add(target_walltime, Scheduler(20), 1.)
        s.run()
        self.assertTrue(s.wall_time() > 1.)

    def test_target_restart_fake(self):
        f = '/tmp/test_simulation/restart/trajectory'
        s = Simulation(DryRun(), output_path=f)
        #s.add(WriterThermo(), Scheduler(20))
        s.add(write_thermo, Scheduler(20))
        s.run(100)
        s.run(100)
        data = numpy.loadtxt(f + '.thermo', unpack=True)
        self.assertEqual(int(data[0][-1]), 200)

    def test_scheduler(self):
        class Simulation:
            def __init__(self):
                self.current_step = 0
        s = Scheduler(3)
        sim = Simulation()
        inext = []
        for i in range(8):
            sim.current_step = i
            inext.append(s(sim))

        self.assertEqual(inext, [3, 3, 3, 6, 6, 6, 9, 9])

    def test_system(self):
        """
        Test that system in Simulation tracks the one in the backend even
        when the latter is reassigned.
        """
        s = Simulation(DryRun(), output_path=None, steps=10)
        s.run()
        s.backend.system = None
        self.assertTrue(s.system is s.backend.system)

    def test_composite(self):
        """
        Test that composite simulation instances (a simulation within a
        simulation object) run independent of their parent instance.

        This checks that there are no regression against the bug fixed
        in 63a7e7863.
        """
        class NewSimulation(Simulation):

            def __init__(self, sim, steps=0, output_path=None, restart=False):
                Simulation.__init__(self, DryRun(), output_path=output_path,
                                    steps=steps, restart=restart)
                self.sim = sim

            def __str__(self):
                return 'NewSimulation'

            def run_until(self, steps):
                self.sim.run()
                self.current_step = steps

        sim = Simulation(DryRun(), steps=3)
        new_sim = NewSimulation(sim, steps=1)
        new_sim.run()
        self.assertEqual(new_sim.current_step, 1)
        self.assertEqual(sim.current_step, 3)

    def test_shell_stop(self):
        from atooms.simulation import shell_stop
        f = '/tmp/test_simulation/shell/trajectory'
        s = Simulation(DryRun(), output_path=f)
        # TODO: why not working?
        #s.add(shell_stop, Scheduler(steps=[20]), 'exit 0')
        s.add(shell_stop, Scheduler(20), 'exit 1')
        s.add(write_thermo, Scheduler(10))
        s.run(100)
        self.assertEqual(s.current_step, 20)

        s = Simulation(DryRun(), output_path=f)
        s.add(shell_stop, Scheduler(20), 'exit 0')
        s.add(write_thermo, Scheduler(10))
        s.run(100)
        self.assertEqual(s.current_step, 100)

        # Test formatted string
        s = Simulation(DryRun(), output_path=f)
        s.add(shell_stop, Scheduler(20), '[ {sim.current_step} -eq 40 ] && exit 1 || exit 0')
        s.run(100)
        self.assertEqual(s.current_step, 40)

    def test_python_stop(self):
        from atooms.simulation import target_python_stop
        f = '/tmp/test_simulation/python/trajectory'
        s = Simulation(DryRun(), output_path=f)
        s.add(target_python_stop, Scheduler(20), '{current_step} == 40')
        s.add(write_thermo, Scheduler(10))
        s.run(100)
        self.assertEqual(s.current_step, 40)

    def tearDown(self):
        rmd('/tmp/test_simulation')


if __name__ == '__main__':
    unittest.main()
