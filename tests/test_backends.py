#!/usr/bin/env python

import unittest
import sys
import os
import numpy

from atooms.simulation import Simulation, target
try:
    import rumd
    from rumdSimulation import rumdSimulation
    SKIP = False
except ImportError:
    SKIP = True

if not SKIP:
    from atooms.backends.rumd import System, Trajectory
    from atooms.backends.rumd import RUMD as Backend

from atooms.core.utils import setup_logging


setup_logging(level=40)

xyz = """\
     3
ioformat=1 dt=0.005000000 boxLengths=6.34960421,6.34960421,6.34960421 numTypes=1 Nose-Hoover-Ps=-0.027281716 Barostat-Pv=0.000000000 mass=1.0000000 columns=type,x,y,z,vx,vy,vz
0       -3.159757      3.145206     -3.145651 1.0 0.0 -1.0
0       -2.986284      3.045374     -2.362755 0.0 1.0 0.0
0       -2.813011      2.380848     -1.037014 -1.0 -1.0 1.0
"""

xyz_2 = """\
     4
ioformat=1 dt=0.005000000 boxLengths=6.34960421,6.34960421,6.34960421 numTypes=2 Nose-Hoover-Ps=-0.027281716 Barostat-Pv=0.000000000 mass=1.0000000,2.000000000 columns=type,x,y,z,vx,vy,vz
0       -3.159757      3.145206     -3.145651 1.0 0.0 1.0
0       -2.986284      3.045374     -2.362755 0.0 1.0 1.0
0       -2.813011      2.380848     -1.037014 0.0 1.0 1.0
1       -1.813011      1.380848     -0.037014 1.0 1.0 0.0
"""

xyz_io2 = """\
4
ioformat=2 timeStepIndex=1234 numTypes=4 integrator=IntegratorNVT,0.00400000019,0.346464646,0.200000003,0 sim_box=RectangularSimulationBox,5.65945244,5.65945244,5.65945244 mass=1.0,0.57,1e20,1e20 columns=type,x,y,z,imx,imy,imz,vx,vy,vz,fx,fy,fz,pe,vir
0 1.859081030 0.200642005 -0.850564003 0 0 0 -0.956519008 -0.441635996 -0.600809991 0.0 0.0 0.0 0.0 0.0
0 1.673092008 0.005931000 0.235958993 0 0 0 -0.661297977 -0.039301001 0.306021005 0.0 0.0 0.0 0.0 0.0
0 -0.317158997 2.685316086 2.080347061 0 0 0 -0.229003996 -0.226451993 0.536197007 0.0 0.0 0.0 0.0 0.0
0 -0.468039989 -0.332848012 1.445088983 0 0 0 -0.344832987 -0.324809015 -0.497487992 0.0 0.0 0.0 0.0 0.0
"""

# TODO: introduce generic tests for backends


class TestBackendRUMD(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('no rumd')

        self.dout = '/tmp/test_backends_rumd_out'
        self.finp = '/tmp/test_backends_rumd_in.xyz'
        with open(self.finp, 'w') as fh:
            fh.write(xyz)

        self.backend = Backend(self.finp)
        self.sim = self.backend.rumd_simulation
        self.backend.rumd_simulation.sample.SetOutputDirectory(self.dout)
        self.backend.rumd_simulation.suppressAllOutput = True
        p = rumd.Pot_LJ_12_6(cutoff_method=rumd.ShiftedPotential)
        p.SetVerbose(False)
        p.SetParams(0, 0, 1., 1., 2.5)
        self.backend.rumd_simulation.SetPotential(p)
        itg = rumd.IntegratorNVT(targetTemperature=2.0, timeStep=0.002)
        self.backend.rumd_simulation.SetIntegrator(itg)

        self.finp2 = '/tmp/test_backends_rumd_in2.xyz'
        with open(self.finp2, 'w') as fh:
            fh.write(xyz_2)
        self.sim2 = rumdSimulation(self.finp2, verbose=False)

        self.finp_io2 = '/tmp/test_backends_rumd_io2.xyz'
        with open(self.finp_io2, 'w') as fh:
            fh.write(xyz_io2)

        from atooms.core.utils import mkdir
        mkdir('/tmp/test_backends')
        self.finp_io2_base = '/tmp/test_backends/0000001.xyz'
        with open(self.finp_io2_base, 'w') as fh:
            fh.write(xyz_io2)

    def test_system(self):
        system = System(self.sim.sample)
        U = system.potential_energy()
        T = system.temperature
        Uref = 36.9236726612
        Tref = 6.0 / (9 - 3)
        # Note places is the number of decimal places, not significant digits, 4 is enough
        self.assertAlmostEqual(U, Uref, 4)
        self.assertAlmostEqual(T, Tref)

    def test_temperature_mass(self):
        system = System(self.sim2.sample)
        T = system.temperature
        Tref = 10.0 / (12 - 3)  # if we don't have the right masses this will fail
        self.assertAlmostEqual(T, Tref)

    def test_particle(self):
        system = System(self.sim.sample)
        p = system.particle
        ref = numpy.array([-3.1597569, 3.14520597, -3.1456511])
        self.assertLess(max(abs(p[0].position - ref)), 1e-6)

    def test_particle_mass(self):
        system = System(self.sim2.sample)
        p = system.particle
        for mref, m in zip(numpy.array([1., 1., 1., 2.]), [pi.mass for pi in p]):
            self.assertAlmostEqual(m, mref)

    def test_trajectory_one_step(self):
        t = Trajectory(self.dout + '_one_step', 'w')
        system = System(self.sim.sample)
        t.write(system, 0)
        t.close()

    def test_trajectory_read(self):
        from atooms.trajectory import TrajectoryRUMD
        with TrajectoryRUMD(self.finp_io2, 'r') as t:
            self.assertEqual(t.steps, [1234])

    def test_trajectory_read_base(self):
        from atooms.trajectory import TrajectoryRUMD
        with TrajectoryRUMD(self.finp_io2_base, 'r') as t:
            self.assertEqual(t.steps, [1])

    def test_simulation(self):
        s = Simulation(self.backend, self.dout, steps=1)
        system = System(self.sim.sample)
        s.system = system
        s.run()

    def test_rmsd(self):
        s = Simulation(self.backend, self.dout, steps=1)
        s.run()
        self.assertGreater(s.rmsd, 0.0)

    def test_target_rmsd(self):
        s = Simulation(self.backend, self.dout, steps=sys.maxsize)
        s.add(target, 10, 'rmsd', 0.3)
        s.run()
        self.assertGreater(s.current_step, 1)
        self.assertGreater(s.rmsd, 0.3)

    def test_checkpoint(self):
        s = Simulation(self.backend, self.dout, steps=10, checkpoint_interval=10)
        s.run()
        # TODO: this will fail, change test for existence of chk file
        # self.assertTrue(os.path.exists(s.trajectory.filename + '.chk'))

    def tearDown(self):
        from atooms.core.utils import rmf, rmd
        rmf('/tmp/test_backends*')
        rmd('/tmp/test_backends')


if __name__ == '__main__':
    unittest.main()
