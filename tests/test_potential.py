import unittest

try:
    import h5py
    from atooms.trajectory import TrajectoryHDF5
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

from atooms.system import System
from atooms.system.cell import Cell
from atooms.interaction import Interaction, PairPotential, CutOff


class PairPotentialTest(unittest.TestCase):

    def test_potential_cs(self):
        def u0ref(rsq):
            return 4*(1/rsq**6 - 1/rsq**3) - 4*(1/2.5**12 - 1/2.5**6)
        def u1ref(rsq):
            return 24*(2/rsq**6 - 1/rsq**3) / rsq
        p = PairPotential('lennard_jones', {'epsilon':1.0, 'sigma':1.0}, [1,1], CutOff('CS', 2.5))
        p.npoints = 10
        rsq, u0, u1 = p.tabulate()
        for r, x, y in zip(rsq[1:], u0[1:], u1[1:]):
            self.assertAlmostEqual(x, u0ref(r))
            self.assertAlmostEqual(y, u1ref(r))

    def test_interacting_system(self):
        p = PairPotential('lennard_jones', {'epsilon':1.0, 'sigma':1.0}, [1,1], CutOff('CS', 2.5))
        i = Interaction(p, 'atomic')
        s = System()
        s.interaction = i

    @unittest.skipIf(not HAS_HDF5, 'no h5py module')
    def test_write_initial_state(self):
        p = [PairPotential('lennard_jones', {'epsilon':1.0, 'sigma':1.0}, [1,1], CutOff('CS', 2.5))]
        i = [Interaction(p, 'atomic')]
        s = System()
        s.cell = Cell([1.0, 1.0, 1.0])
        t = TrajectoryHDF5('/tmp/test.h5', 'w')
        t.write_interaction(i)
        t.close()

        t = TrajectoryHDF5('/tmp/test.h5', 'r')
        i = t.read_interaction()
        t.close()

if __name__ == '__main__':
    unittest.main()

        
