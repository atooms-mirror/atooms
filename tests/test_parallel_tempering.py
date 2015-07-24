#!/usr/bin/env python

import os
import unittest

from atooms.simulation.parallel_tempering import ParallelTempering

class TestReplicaExchange(unittest.TestCase):

    def test_rx_atooms(self):
        # TODO: refactor replica exchange as a factory, of which the actual adapters define the backend. The only problem is that the backend always require some fine tuning.
        try:
            from atooms.adapters.atoomsf90 import SimulationAtooms
            fref = 'reference/kalj-small.h5'
            if not os.path.exists(fref):
                raise
        except ImportError:
            self.skipTest('no atooms')
        except:
            self.skipTest('no ref')

        seed = 10
        nsteps = 3000
        swap_period = 1000
        Tl = [1.0, 1.1, 1.2]
        dir_root='/tmp'
        file_log = dir_root + '/rx.log'
        file_state_config = [dir_root + '/rx.s%d.h5' % i for i in range(len(Tl))]
        file_replica_out = [dir_root + '/rx.r%d.out' % i for i in range(len(Tl))]

        sim_adapter = [SimulationAtooms(f, file_input=fref,
                                        opts={'--thermostat':'Berendsen',
                                              '--thermostat-temperature':T,
                                              '--dt':0.004}) 
                       for f, T in zip(file_state_config, Tl)]
        rx = ReplicaExchange(Tl, file_state_config, file_log, file_replica_out)
        rx_sim = ReplicaExchangeSimulation(rx, sim_adapter, swap_period)
        rx_sim.run(nsteps, 100)

        ref_all = [[0, 1, 2, 2], [1, 0, 0, 1], [2, 2, 1, 0]]
        for f, ref in zip(file_replica_out, ref_all):
            with open(f) as fh:
                r = [int(l.split()[2]) for l in fh.readlines()]
            self.assertEqual(r, ref)

    def test_rx_rumd(self):

        try:
            import rumd
            from rumdSimulation import rumdSimulation
            from atooms.adapters.rumd import Simulation, System, Trajectory
        except ImportError:
            self.skipTest('no rumd')

        def potential():
            pot = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
            pot.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
            pot.SetParams(i=1, j=0, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
            pot.SetParams(i=0, j=1, Epsilon=1.5, Sigma=0.8, Rcut=2.5)
            pot.SetParams(i=1, j=1, Epsilon=0.5, Sigma=0.88, Rcut=2.5)
            return pot

        def rx(Tl, dir_root, input_file, nsteps, swap_period):
            
            # Create simulation and integrators
            sim = [rumdSimulation(f) for f in input_file]
            igt = [rumd.IntegratorNVT(targetTemperature=T, timeStep=0.002) for T in Tl]
            for s, i in zip(sim, igt):
                s.AddPotential(potential())
                s.SetIntegrator(i)

            dir_state_out = [dir_root + '/r%d' % i for i in range(len(Tl))]
            sim_adapter = [Simulation(s, d) for s, d in zip(sim, dir_state_out)]
            rx_sim = ParallelTempering(dir_root, dir_state_out, Tl, sim_adapter, swap_period)
            rx_sim.setup(target_steps = nsteps,
                         thermo_period = 1,
                         config_period = 1)
            rx_sim.run()
                 
        Tl = [1.0, 1.1, 1.2]
        input_file = ['reference/kalj256.h5.xyz' for T in Tl]
        dir_root = '/tmp/rx'
        rx(Tl, dir_root, input_file, nsteps=3, swap_period=1000)


if __name__ == '__main__':
    unittest.main()


