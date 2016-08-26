import unittest

try:
    import h5py
    from atooms.trajectory import TrajectoryHDF5
    HAS_HDF5 = True
except:
    HAS_HDF5 = False

from atooms.system import System
from atooms.system.cell import Cell
from atooms.system.particle import Particle
from atooms.interaction.interaction import Interaction
from atooms.potential.potential import PairPotential
from atooms.potential.cutoff import CutOff

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
        s.particle = [Particle(position=[1.0, 1.0, 1.0], velocity=[0.0, 0.0, 0.0])]
        s.cell = Cell([1.0, 1.0, 1.0])
        with TrajectoryHDF5('/tmp/test.h5', 'w') as t:
            t.write_interaction(i)
            t.write(s, 0)

        with TrajectoryHDF5('/tmp/test.h5', 'r') as t:
            i = t.read_interaction()
            s = t[0]

    @unittest.skipIf(not HAS_HDF5, 'no h5py module')
    def test_fmt(self):
        with TrajectoryHDF5('/tmp/test.h5', 'w') as t:
            t.exclude(['velocity'])
            t.include(['position'])
        

if __name__ == '__main__':
    unittest.main()

        
