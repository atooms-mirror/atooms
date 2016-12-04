import unittest

try:
    import h5py
    from atooms.trajectory import TrajectoryHDF5
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

from atooms.system import System
from atooms.system.cell import Cell
from atooms.interaction.interaction import Interaction
from atooms.potential.lennard_jones import LennardJones
from atooms.potential.potential import PairPotential
from atooms.potential.cutoff import CutOff

try:
    import h5py
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

class PairPotentialTest(unittest.TestCase):

    def test_potential_cs(self):
        p = LennardJones('LJ', {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))
        p.npoints = 10
        rsq, u0, u1 = p.tabulate()
        self.assertAlmostEqual(u0[2], -0.870024506431)
        self.assertAlmostEqual(u1[2], -1.85584204277)

    def test_interacting_system(self):
        p = LennardJones('LJ', {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))
        i = Interaction("atomic", p)
        s = System()
        s.interaction = i

    @unittest.skipIf(not HAS_HDF5, 'no h5py module')
    def test_write_initial_state(self):
        p = [LennardJones('LJ', {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))]
        i = [Interaction("atomic", p)]
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

        
