#!/usr/bin/env python

import os
import unittest
from atooms.trajectory import TrajectoryXYZ
from atooms.simulation import Simulation, write_thermo, write_config, target
from atooms.core.utils import setup_logging, mkdir, rmd, rmf
try:
    from atooms.backends.lammps import LAMMPS, Interaction
    SKIP = False
except ImportError:
    SKIP = True

setup_logging(level=40)


class Test(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('missing LAMMPS')
        self.input_file = os.path.join(os.path.dirname(__file__),
                                       '../data/lj_N1000_rho1.0.xyz')

    def test_single(self):
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
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.62635, places=5)
        sim.run(10)
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.64526, places=5)
        sim.run(10)
        x = sim.system.particle[0].position[0]
        self.assertAlmostEqual(x, 3.675987, places=5)

    def test_nvt(self):
        import sys
        import random
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nvt temp 2.0 2.0 0.2
        timestep        0.002
        """
        def store(sim, T, U):
            T.append(sim.system.temperature)
            U.append(sim.system.potential_energy(per_particle=True))

        with TrajectoryXYZ(self.input_file) as th:
            system = th[-1]

        for inp in [self.input_file, TrajectoryXYZ(self.input_file),
                    system]:
            T, U = [], []
            random.seed(1)
            bck = LAMMPS(inp, cmd)
            sim = Simulation(bck)
            sim.system.temperature = 1.5
            sim.add(store, 500, T, U)
            sim.run(2000)
            ave = sum(T[3:]) / len(T[3:])
            self.assertAlmostEqual(ave, 2.0, places=1)
            if isinstance(inp, TrajectoryXYZ):
                inp.close()

    def test_nvt_nofix(self):
        import sys
        import random
        from atooms.system import Thermostat
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        timestep        0.002
        """
        random.seed(1)
        T = []
        def store(sim, T):
            T.append(sim.system.temperature)
        bck = LAMMPS(self.input_file, cmd)
        sim = Simulation(bck)
        sim.system.temperature = 1.4
        sim.system.thermostat = Thermostat(2.0)
        sim.add(store, 500, T)
        sim.run(4000)
        ave = sum(T[3:]) / len(T[3:])
        self.assertAlmostEqual(ave, 2.0, places=1)

    def test_energy(self):
        import sys
        import numpy
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        """
        bck = LAMMPS(self.input_file, cmd)
        bck.system.interaction.compute("energy", bck.system.particle, bck.system.cell)
        self.assertEqual(bck.system.interaction.energy / len(bck.system.particle), -4.2446836)  # crosschecked

        # Relaxed FCC
        bck = LAMMPS(os.path.join(os.path.dirname(__file__),
                                  '../data/lj_fcc_N108.xyz'), cmd)
        bck.system.interaction.compute("forces", bck.system.particle, bck.system.cell)
        # Test norm of force per particle
        U = bck.system.potential_energy(per_particle=True)
        W = bck.system.force_norm(per_particle=True)
        P = bck.system.pressure
        self.assertAlmostEqual(U, -7.7615881, places=7)
        self.assertAlmostEqual(W, 4.2e-11, places=2)
        self.assertAlmostEqual(P, -3.3935748, places=7)

    def test_trajectory(self):
        import sys
        with open('/tmp/test_lammps.atom', 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-3 3
-3 3
-3 3
ITEM: ATOMS id type xs ys zs
2 1 0.10 0.11 0.12
1 1 0.20 0.21 0.22
ITEM: TIMESTEP
1
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-4 4
-4 4
-4 4
ITEM: ATOMS id type xs ys zs
1 1 0.00 0.01 0.02
2 1 0.50 0.51 0.52
""")
        from atooms.trajectory import TrajectoryLAMMPS
        def scale(pos, side):
            return [(x - 0.5) * L for x, L in zip(pos, side)]
        with TrajectoryLAMMPS('/tmp/test_lammps.atom') as th:
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
            self.assertEqual(list(th[0].particle[0].position), scale([0.20, 0.21, 0.22], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[0].particle[1].position), scale([0.10, 0.11, 0.12], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[1].cell.side), [8.0, 8.0, 8.0])
            self.assertEqual(list(th[1].particle[0].position), scale([0.00, 0.01, 0.02], [8.0, 8.0, 8.0]))
            self.assertEqual(list(th[1].particle[1].position), scale([0.50, 0.51, 0.52], [8.0, 8.0, 8.0]))

    def test_trajectory_folder(self):
        import sys
        mkdir('/tmp/test_lammps.d')
        with open('/tmp/test_lammps.d/0.atom', 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-3 3
-3 3
-3 3
ITEM: ATOMS id type xs ys zs
2 1 0.10 0.11 0.12
1 1 0.20 0.21 0.22
""")
        with open('/tmp/test_lammps.d/1.atom', 'w') as fh:
            fh.write("""\
ITEM: TIMESTEP
20
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
-4 4
-4 4
-4 4
ITEM: ATOMS id type xs ys zs
1 1 0.00 0.01 0.02
2 1 0.50 0.51 0.52
""")
        from atooms.trajectory import TrajectoryFolderLAMMPS
        def scale(pos, side):
            return [(x - 0.5) * L for x, L in zip(pos, side)]
        with TrajectoryFolderLAMMPS('/tmp/test_lammps.d') as th:
            self.assertEqual(th.steps, [10, 20])
            self.assertEqual(list(th[0].cell.side), [6.0, 6.0, 6.0])
            self.assertEqual(list(th[0].particle[0].position), scale([0.20, 0.21, 0.22], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[0].particle[1].position), scale([0.10, 0.11, 0.12], [6.0, 6.0, 6.0]))
            self.assertEqual(list(th[1].cell.side), [8.0, 8.0, 8.0])
            self.assertEqual(list(th[1].particle[0].position), scale([0.00, 0.01, 0.02], [8.0, 8.0, 8.0]))
            self.assertEqual(list(th[1].particle[1].position), scale([0.50, 0.51, 0.52], [8.0, 8.0, 8.0]))

    # def test_species(self):
    #     import sys
    #     import atooms.trajectory
    #     cmd = """
    #     pair_style      lj/cut 2.5
    #     pair_coeff      1 1 1.0 1.0 2.5
    #     fix             1 all nve
    #     """
    #     bck = LAMMPS(os.path.join(os.path.dirname(__file__),
    #                               '../data/trajectory.xyz'), cmd)
    #     sim = Simulation(bck)
    #     xi = sim.system.particle[0].species
    #     sim.run(1)
    #     xf = sim.system.particle[0].species
    #     print xi, xf

    def test_roundoff(self):
        import sys
        cmd = """
        pair_style      lj/cut 2.5
        pair_coeff      1 1 1.0 1.0 2.5
        neighbor        0.3 bin
        neigh_modify    every 20 delay 0 check no
        fix             1 all nve
        """
        import numpy
        numpy.set_printoptions(precision=15)
        bck = LAMMPS(self.input_file, cmd)
        bck.run(10)
        x1 = bck.system.particle[3].position[0]
        bck.run(0)
        x2 = bck.system.particle[3].position[0]
        # The difference should be of order of machine precision
        # Unfortunately, without packing/unpacking data in binary
        # format, we cannot maintain coehrence at machine precision
        self.assertAlmostEqual(abs(x2-x1), 0.0, places=12)

    def tearDown(self):
        rmd('/tmp/test_lammps.d')
        rmf('/tmp/test_lammps*')


if __name__ == '__main__':
    unittest.main()
