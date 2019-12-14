#!/usr/bin/env python

import unittest

import numpy
from atooms.system import System, Particle, Cell
from atooms.trajectory.ram import TrajectoryRam, TrajectoryRamFull

class Test(unittest.TestCase):

    def setUp(self):
        pass

    @unittest.expectedFailure
    def test_ram_inplace(self):
        particle = [Particle([0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRam()
        t[0] = system
        particle[0].position += 1.0
        self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram(self):
        particle = [Particle([0.0, 0.0, 0.0])]
        system = System(particle)
        t = TrajectoryRam()
        t[0] = system
        particle[0].position = numpy.array([1.0, 1.0, 1.0])
        self.assertFalse((t[0].particle[0].position == particle[0].position).all())

    def test_ram_full(self):
        particle = [Particle([0.0, 0.0, 0.0])]
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
        from atooms.backends.lammps import LAMMPS
        from atooms.simulation import Simulation
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/lj_N1000_rho1.0.xyz')
        import sys
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        bck = LAMMPS(self.input_file, cmd)
        sim = Simulation(bck)
        sim.system.particle[0].position = numpy.zeros(3)

        t = TrajectoryRamFull()
        t[0] = sim.system
        sim.system.particle[0].position += 1.0
        self.assertFalse((sim.system.particle[0].position == t[0].particle[0].position).all())

    def test_ram_unfolded(self):
        import os
        import numpy
        from atooms.backends.lammps import LAMMPS
        from atooms.trajectory import TrajectoryXYZ, Unfolded
        input_file = os.path.join(os.path.dirname(__file__),
                                  '../data/lj_N1000_rho1.0.xyz')
        tinp = TrajectoryXYZ(input_file)
        tf = TrajectoryRamFull()
        t = TrajectoryRam()
        for i, s in enumerate(tinp):
            tf[i] = s
            t[i] = s
        self.assertTrue((tf[-1].particle[-1].position == tinp[-1].particle[-1].position).all())
        self.assertTrue((t[-1].particle[-1].position == tinp[-1].particle[-1].position).all())

        tu = Unfolded(t, fixed_cm=False)
        for s in tu:
            pass
        self.assertTrue((t[-1].particle[-1].position == tinp[-1].particle[-1].position).all())

        tu = Unfolded(tf, fixed_cm=False)
        for s in tu:
            pass
        self.assertTrue((tf[-1].particle[-1].position == tinp[-1].particle[-1].position).all())

        tu = Unfolded(t, fixed_cm=True)
        for s in tu:
            pass
        self.assertTrue((t[-1].particle[-1].position == tinp[-1].particle[-1].position).all())

        tu = Unfolded(tf, fixed_cm=True)
        for s in tu:
            pass
        self.assertTrue((tf[-1].particle[-1].position == tinp[-1].particle[-1].position).all())


if __name__ == '__main__':
    unittest.main()


