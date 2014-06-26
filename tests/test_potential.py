import unittest

from atooms.system import System
from atooms.system.cell import Cell
from atooms.trajectory import TrajectoryHDF5
from atooms.interaction.interaction import Interaction
from atooms.potential.potential import PairPotential
from atooms.potential.cutoff import CutOff

try:
    import h5py
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

class PairPotentialTest(unittest.TestCase):

    def test_potential(self):
        p = PairPotential("LennardJones", {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))

    def test_interacting_system(self):
        p = [PairPotential("LennardJones", {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))]
        i = Interaction("atomic", p)
        s = System()
        s.interaction = i

    @unittest.skipIf(not HAS_HDF5, 'no h5py module')
    def test_write_initial_state(self):
        p = [PairPotential("LennardJones", {"epsilon":1.0, "sigma":1.0}, [1,1], CutOff("CS", 2.5))]
        i = [Interaction("atomic", p)]
        s = System()
        s.cell = Cell()
        #s.interaction = i
        t = TrajectoryHDF5('/tmp/test.h5', 'w')
        #t.write_initial_state(s)
        t.write_initial_state_interaction(i)
        t.close()

        t = TrajectoryHDF5('/tmp/test.h5', 'r')
        #t.write_initial_state(s)
        i = t.read_initial_state_interaction()
        t.close()

        #print i[0].potential[0]

if __name__ == '__main__':
    unittest.main()

        
