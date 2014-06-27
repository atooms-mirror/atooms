#!/usr/bin/env python

import unittest
import os
import numpy

import rumd
from rumdSimulation import rumdSimulation
from atooms.adapters.rumd import Simulation, System, Trajectory

xyz = """\
     4
ioformat=1 dt=0.005000000 boxLengths=6.34960421,6.34960421,6.34960421 numTypes=2 Nose-Hoover-Ps=-0.027281716 Barostat-Pv=0.000000000 mass=1.0000000,2.000000000 columns=type,x,y,z,vx,vy,vz
0       -3.159757      3.145206     -3.145651 1.0 0.0 1.0
0       -2.986284      3.045374     -2.362755 0.0 1.0 1.0
0       -2.813011      2.380848     -1.037014 0.0 1.0 1.0
1       -2.813011      2.380848     -1.037014 1.0 1.0 0.0
"""

# TODO: make test_adapters a package

class TestAdaptersRUMD(unittest.TestCase):

    def setUp(self):
        self.fout = '/tmp/test_adapter_rumd_out.xyz.gz'
        self.dout = '/tmp/test_adapter_rumd_out'
        self.finp = '/tmp/test_adapter_rumd_in.xyz'
        fh = open(self.finp, 'w')
        fh.write(xyz)
        fh.close()

        self.s = rumdSimulation(self.finp, verbose=False)
        self.s.SetOutputScheduling("energies", "none")
        self.s.SetOutputScheduling("trajectory", "none")
        self.s.sample.SetOutputDirectory(self.dout)
        self.s.suppressAllOutput = True
        p = rumd.Pot_LJ_12_6(cutoff_method = rumd.ShiftedPotential)
        p.SetVerbose(False)
        p.SetParams(0, 0, 1., 1., 2.5)
        self.s.SetPotential(p)
        i = rumd.IntegratorNVT(targetTemperature=2.0, timeStep=0.002)
        self.s.SetIntegrator(i)

    def test_system(self):
        system = System(self.s)
        U = system.potential_energy()
        T = system.temperature()
        Uref = 36.9236726612
        Tref = 20.0/9
        self.assertLess(abs((U-Uref)/Uref), 0.001)
        self.assertLess(abs((T-Tref)/Tref), 1e-9)

    def test_particle(self):
        system = System(self.s)
        p = system.particle
        ref = numpy.array([-3.1597569, 3.14520597, -3.1456511])
        mref = numpy.array([1.,1.,1.,2.])
        self.assertLess(max(abs(p[0].position - ref)), 1e-6)
        self.assertLess(max(abs([pi.mass for pi in p] - mref)), 1e-6)

    def test_trajectory(self):
        t = Trajectory(self.fout, 'w')
        system = System(self.s)
        t.write_sample(system, 0, 0)
        t.close()

    def test_simulation(self):
        s = Simulation(self.s, self.dout)
        s.setup(target_steps = 1)
        system = System(self.s)
        s.system = system
        s.run()

    def test_rmsd(self):
        s = Simulation(self.s, self.dout)
        s.setup(target_steps = 1)
        s.run()
        self.assertGreater(s.rmsd, 0.0)

    def test_target_rmsd(self):
        s = Simulation(self.s, self.dout)       
        s.setup(target_rmsd = 0.3)
        s.run()
        self.assertGreater(s.steps, 1)
        self.assertGreater(s.rmsd, 0.3)

    def test_checkpoint(self):
        s = Simulation(self.s, self.dout)
        s.setup(target_steps = 10, checkpoint_period = 10)
        s.run()
        self.assertTrue(os.path.exists(s.trajectory.filename + '.chk'))
        s.setup(target_steps = 20, checkpoint_period = 10, reset=True)
        s.restart = True
        s.run()
                
    def tearDown(self):
        import shutil
        if os.path.exists(self.finp):
            os.remove(self.finp)
        if os.path.exists(self.fout):
            os.remove(self.fout)
        if os.path.exists(self.dout):
            shutil.rmtree(self.dout)

if __name__ == '__main__':
    # import logging
    # l = logging.getLogger()
    # l.setLevel('DEBUG')
    unittest.main(verbosity=0)


