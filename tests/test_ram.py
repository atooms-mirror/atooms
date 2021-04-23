#!/usr/bin/env python

import unittest

import numpy
from atooms.system import System, Particle, Cell
from atooms.trajectory.ram import TrajectoryRam, TrajectoryRamFull

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_ram_inplace(self):
        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRam()
        t[0] = system
        particle[0].position += 1.0
        self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram(self):
        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRam()
        t[0] = system
        particle[0].position = numpy.array([1.0, 1.0, 1.0])
        self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram_full(self):
        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRamFull()
        t[0] = system
        system.particle[0].position = numpy.array([2.0, 2.0, 2.0])
        t[0] = system
        particle[0].position += 1.0
        self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram_lammps(self):
        import os
        import numpy
        import sys
        from atooms.backends.lammps import LAMMPS
        from atooms.simulation import Simulation
        import atooms.backends.lammps
        if not atooms.backends.lammps.installed():
            self.skipTest('missing LAMMPS')
        input_file = os.path.join(os.path.dirname(__file__),
                                  '../data/lj_N1000_rho1.0.xyz')
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        bck = LAMMPS(input_file, cmd)
        sim = Simulation(bck)
        sim.system.particle[0].position = numpy.zeros(3)

        t = TrajectoryRamFull()
        t[0] = sim.system
        sim.system.particle[0].position += 1.0
        self.assertFalse((sim.system.particle[0].position == t[0].particle[0].position).all())

    def test_ram_unfolded(self):
        import os
        import numpy
        import random
        from atooms.backends.lammps import LAMMPS
        from atooms.trajectory import TrajectoryXYZ, Unfolded
        from atooms.simulation import Simulation
        from atooms.system import Thermostat
        import atooms.backends.lammps
        if not atooms.backends.lammps.installed():
            self.skipTest('missing LAMMPS')

        def store(sim, ram):
            ram.write(sim.system, step=sim.current_step)

        tf = TrajectoryRamFull()
        t = TrajectoryRam()
        with TrajectoryXYZ(os.path.join(os.path.dirname(__file__),
                                        '../data/ka_N150_rho1.200.xyz')) as tinp:
            for i, s in enumerate(tinp):
                tf[i] = s
                t[i] = s
        tuf = Unfolded(tf, fixed_cm=True)
        tu = Unfolded(t, fixed_cm=True)
        for s, sf in zip(tu, tuf):
            for p, pf in zip(s.particle, sf.particle):
                self.assertTrue((p.position == pf.position).all())

        def mobility(t):
            import numpy as np
            side = t[0].cell.side
            K = 0
            unfoldedtj = Unfolded(t, fixed_cm=True)
            pos_0 = unfoldedtj[0].dump('pos')
            for j in range(1, len(t)):
                pos_1 = unfoldedtj[j].dump('pos')
                K += np.sum((pos_1 - pos_0)**2)
                pos_0 = pos_1
            return K

        # print mobility(t), mobility(tf)

    def test_ram_copy(self):
        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRam()
        t[0] = system
        t[0].particle[0].position = numpy.array([1.0, 1.0, 1.0])
        # print system.particle[0].position, t[0].particle[0].position
        # print id(t[0].particle[0])
        # print id(t[0].particle[0])

        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        s = System(particle)
        t = TrajectoryRamFull()
        t[0] = system
        s.update(t[0])
        s.particle[0].position = numpy.array([1.0, 1.0, 1.0])
        # print system.particle[0].position, t[0].particle[0].position, s.particle[0].position
        # print id(t[0].particle[0])
        # print id(t[0].particle[0])

        particle = [Particle(position=[0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRamFull()
        t[0] = system
        system.particle[0].position = numpy.array([1.0, 1.0, 1.0])
        # print system.particle[0].position, t[0].particle[0].position
        # print id(t[0].particle[0])
        # print id(system.particle[0])

        #self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram_lammps_write(self):
        import os
        import numpy
        import sys
        from atooms.backends.lammps import LAMMPS
        from atooms.simulation import Simulation
        from atooms.simulation.observers import write_to_ram
        import atooms.backends.lammps
        if not atooms.backends.lammps.installed():
            self.skipTest('missing LAMMPS')

        input_file = os.path.join(os.path.dirname(__file__),
                                  '../data/lj_N1000_rho1.0.xyz')
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        bck = LAMMPS(input_file, cmd)
        ram = TrajectoryRamFull()
        sim = Simulation(bck)
        sim.add(write_to_ram, 10, ram)
        sim.run(100)
        self.assertEqual(len(ram.steps), 11)

    def test_callback_copy(self):
        import copy
        from atooms.trajectory.decorators import filter_species
        particle = [Particle(species='A'), Particle(species='B')]
        system = System(particle)
        t = TrajectoryRamFull()
        t[0] = system
        t.add_callback(copy.deepcopy)
        t.add_callback(filter_species, 'A')
        self.assertEqual(t[0].distinct_species, ['A'])
        t.callbacks.pop()
        t.callbacks.pop()
        self.assertEqual(t[0].distinct_species, ['A', 'B'])

if __name__ == '__main__':
    unittest.main()
