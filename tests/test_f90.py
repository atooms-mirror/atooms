import os
import unittest

try:
    from atooms.backends import f90
    HAS_F90 = True
except ImportError:
    HAS_F90 = False
    
try:
    import atooms.models
    HAS_MODELS = True
except ImportError:
    HAS_MODELS = False
    
from atooms.trajectory import Trajectory


class Test(unittest.TestCase):

    def setUp(self):
        if not HAS_F90:
            self.skipTest("skip f90 backend tests (missing f2py_jit?)")
        self.trajectory = f90.Trajectory('data/lj_N256_rho1.0.xyz')
        #self.trajectory = f90.Trajectory('data/lj_N4000.xyz')

    @unittest.skipIf(not HAS_MODELS, 'no atooms-models module')
    def test_collinear(self):
        import atooms.models
        from atooms.backends.f90 import Interaction, System, Particle, Cell
        db = atooms.models.load()
        model = db["lennard_jones"]
        particles = [Particle(position=[0.0, 0.0, 0.0], species=1),
	             Particle(position=[1.0, 0.0, 0.0], species=1),
	             Particle(position=[2.0, 0.0, 0.0], species=1)]
        cell = Cell([10., 10., 10.])	    
        system = System(particles, cell)
        system.interaction = Interaction(model)
        self.assertAlmostEqual(system.potential_energy(), -0.01257276409199999)
        
    def test_trajectory(self):
        # Make sure the original trajectories have no wrappers
        with Trajectory('data/lj_N256_rho1.0.xyz') as trajectory:
            self.assertEqual(trajectory.class_callbacks, None)

    @unittest.skipIf(not HAS_MODELS, 'no atooms-models module')
    def test_energy(self):
        system = self.trajectory[0]
        system.compute_interaction()
        self.assertAlmostEqual(system.potential_energy(per_particle=True, cache=True), -3.8079776291909284)

    def test_neighbors(self):
        rcut = [[1.5]]
        system = self.trajectory[0]
        system.neighbor_list = f90.NeighborList(rcut)
        system._compute_neighbor_list()
        self.assertEqual(list(system.particle[0].neighbors), [17, 32, 52, 91, 109, 112, 121, 140, 149, 162, 181, 198, 231, 247, 256])

    @unittest.skipIf(not HAS_MODELS, 'no atooms-models module')
    def test_model(self):
        with f90.Trajectory('data/lj_rho0.85.xyz') as trajectory:
            trajectory.metadata['model'] = 'lennard_jones'
            system = trajectory[0]
            self.assertAlmostEqual(system.potential_energy(per_particle=True), -2.24379330538)

    def tearDown(self):
        self.trajectory.close()
        
if __name__ == '__main__':
    unittest.main()
        

