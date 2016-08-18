#!/usr/bin/env python

import unittest
import sys
import os
import numpy
from atooms.simulation import log

log.setLevel(40)

from atooms.simulation import Simulation
try:
    import rumd
    from rumdSimulation import rumdSimulation
    from atooms.backends.rumd_backend import System, Trajectory
    from atooms.backends.rumd_backend import RumdBackend as Backend
    SKIP = False
except ImportError:
    SKIP = True


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

# TODO: make test_backend a package

class TestBackendRUMD(unittest.TestCase):

    def setUp(self):
        if SKIP:
            self.skipTest('no rumd')

        self.fout = '/tmp/test_adapter_rumd_out.xyz.gz'
        self.dout = '/tmp/test_adapter_rumd_out'
        self.finp = '/tmp/test_adapter_rumd_in.xyz'
        with open(self.finp, 'w') as fh:
            fh.write(xyz)

        self.sim = rumdSimulation(self.finp, verbose=False)
        self.sim.SetOutputScheduling("energies", "none")
        self.sim.SetOutputScheduling("trajectory", "none")
        self.sim.sample.SetOutputDirectory(self.dout)
        self.sim.suppressAllOutput = True
        p = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
        p.SetVerbose(False)
        p.SetParams(0, 0, 1., 1., 2.5)
        self.sim.SetPotential(p)
        i = rumd.IntegratorNVT(targetTemperature=2.0, timeStep=0.002)
        self.sim.SetIntegrator(i)

        self.finp2 = '/tmp/test_adapter_rumd_in2.xyz'
        with open(self.finp2, 'w') as fh:
            fh.write(xyz_2)
        self.sim2 = rumdSimulation(self.finp2, verbose=False)

    def test_system(self):
        system = System(self.sim.sample)
        U = system.potential_energy()
        T = system.temperature()
        Uref = 36.9236726612
        Tref = 2*6.0/6
        # Note places is the number of decimal places, not significant digits, 4 is enough
        self.assertAlmostEqual(U, Uref, 4)
        self.assertAlmostEqual(T, Tref)

    def test_temperature_mass(self):
        system = System(self.sim2.sample)
        T = system.temperature()
        Tref = 20.0/9 # if we don't have the right masses this will fail
        self.assertAlmostEqual(T, Tref)

    def test_particle(self):
        system = System(self.sim.sample)
        p = system.particle
        ref = numpy.array([-3.1597569, 3.14520597, -3.1456511])
        self.assertLess(max(abs(p[0].position - ref)), 1e-6)

    def test_particle_mass(self):
        system = System(self.sim2.sample)
        p = system.particle
        for mref, m in zip(numpy.array([1.,1.,1.,2.]), [pi.mass for pi in p]):
            self.assertAlmostEqual(m, mref)

    def test_trajectory(self):
        t = Trajectory(self.fout, 'w')
        system = System(self.sim.sample)
        t.write(system, 0)
        t.close()

    def test_simulation(self):
        s = Simulation(Backend(self.sim), self.dout, steps = 1)
        system = System(self.sim.sample)
        s.system = system
        s.run()

    def test_rmsd(self):
        s = Simulation(Backend(self.sim), self.dout, steps = 1)
        s.run()
        self.assertGreater(s.rmsd, 0.0)

    def test_target_rmsd(self):
        s = Simulation(Backend(self.sim), self.dout, rmsd = 0.3, steps=sys.maxint)
        s.run()
        self.assertGreater(s.steps, 1)
        self.assertGreater(s.rmsd, 0.3)

    def test_checkpoint(self):
        s = Simulation(Backend(self.sim), self.dout, steps = 10, checkpoint_interval = 10)
        s.run()
        # TODO: this will fail, change test for existence of chk file
        # self.assertTrue(os.path.exists(s.trajectory.filename + '.chk'))
        # TODO: this will fail, because changing steps should update scheduler!
        # s.target_steps = 20
        # s.restart = True
        # s.run()
                
    def tearDown(self):
        import shutil
        if os.path.exists(self.finp):
            os.remove(self.finp)
        if os.path.exists(self.finp2):
            os.remove(self.finp2)
        if os.path.exists(self.fout):
            os.remove(self.fout)
        if os.path.exists(self.dout):
            shutil.rmtree(self.dout)

if __name__ == '__main__':
    # import logging
    # l = logging.getLogger()
    # l.setLevel('DEBUG')
    unittest.main(verbosity=0)


