import unittest

try:
    import h5py
    from atooms.trajectory import TrajectoryHDF5
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

from atooms.core.utils import rmf
from atooms.system import System
from atooms.system.cell import Cell
from atooms.interaction import Interaction, PairPotential, CutOff


class PairPotentialTest(unittest.TestCase):

    def test_potential_cs(self):
        def u0ref(rsq):
            return 4 * (1/rsq**6 - 1/rsq**3) - 4 * (1/2.5**12 - 1/2.5**6)

        def u1ref(rsq):
            return 24 * (2/rsq**6 - 1/rsq**3) / rsq
        p = PairPotential('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0},
                          [1, 1], CutOff('CS', 2.5))
        p.npoints = 10
        rsq, u0, u1 = p.tabulate()
        for r, x, y in zip(rsq[1:], u0[1:], u1[1:]):
            self.assertAlmostEqual(x, u0ref(r))
            self.assertAlmostEqual(y, u1ref(r))

    def test_interacting_system(self):
        p = PairPotential('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0},
                          [1, 1], CutOff('CS', 2.5))
        i = Interaction(p, 'atomic')
        s = System()
        s.interaction = i

    @unittest.skipIf(not HAS_HDF5, 'no h5py module')
    def test_write_initial_state(self):
        p = [PairPotential('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0},
                           [1, 1], CutOff('CS', 2.5))]
        i = [Interaction(p, 'atomic')]
        s = System()
        s.cell = Cell([1.0, 1.0, 1.0])
        t = TrajectoryHDF5('/tmp/test_potential.h5', 'w')
        t.write_interaction(i)
        t.close()

        t = TrajectoryHDF5('/tmp/test_potential.h5', 'r')
        i = t.read_interaction()
        t.close()

    def test_hard_sphere(self):
        p = PairPotential('hard_sphere', {'sigma': 1.0}, [1, 1])
        self.assertEqual(p.hard_core, 1.0)
        self.assertEqual(p.cutoff.radius, 0.0)

    def test_square_well(self):
        p = PairPotential('square_well', {'epsilon': -1.0, 'delta': 0.1, 'sigma': 1.0}, [1, 1])
        self.assertEqual(p.params['epsilon'], -1.0)
        self.assertEqual(p.hard_core, 1.0)
        self.assertEqual(p.cutoff.radius, 1.1)
        self.assertEqual(p.compute(1.05)[0], -1.0)
        self.assertTrue(p.is_zero(1.2**2))

    def test_tabulate_fmt(self):
        from atooms.interaction.potential import tabulate
        txt = tabulate('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0}, npoints=6)
        txt = tabulate('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0}, fmt='uwh', npoints=6)
        txt = tabulate('lennard_jones', {'epsilon': 1.0, 'sigma': 1.0}, fmt='other', npoints=6)
        
    def tearDown(self):
        rmf('/tmp/test_potential.h5')


if __name__ == '__main__':
    unittest.main()
