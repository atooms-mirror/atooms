#!/usr/bin/env python

import random
import copy
import numpy
import unittest
from atooms.system import System
from atooms.system.cell import Cell
from atooms.system.particle import Particle
from atooms.system.reservoir import Thermostat, Barostat, Reservoir


class Test(unittest.TestCase):

    def setUp(self):
        N = 100
        L = 10.0
        self.ref = System()
        self.ref.cell = Cell([L, L, L])
        self.ref.particle = []
        self.ref.thermostat = Thermostat(1.0)
        self.ref.barostat = Barostat(1.0)
        self.ref.reservoir = Reservoir(1.0)
        while len(self.ref.particle) <= N:
            pos = [(random.random() - 0.5) * L,
                   (random.random() - 0.5) * L,
                   (random.random() - 0.5) * L]
            self.ref.particle.append(Particle(position=pos))

    def test_density(self):
        system = copy.copy(self.ref)
        density_old = system.density
        system.density = density_old * 1.1
        self.assertAlmostEqual(system.density, density_old * 1.1)

    def test_temperature(self):
        system = copy.copy(self.ref)
        system.set_temperature(1.0)
        self.assertAlmostEqual(system.temperature, 1.0)

    def test_cm(self):
        system = copy.copy(self.ref)
        system.set_temperature(1.0)
        system.fix_momentum()
        self.assertAlmostEqual(system.cm_velocity[0], 0.0)
        self.assertAlmostEqual(system.cm_velocity[1], 0.0)
        self.assertAlmostEqual(system.cm_velocity[2], 0.0)
        pos_cm = system.cm_position
        for p in system.particle:
            p.position -= pos_cm
        self.assertAlmostEqual(system.cm_position[0], 0.0)
        self.assertAlmostEqual(system.cm_position[1], 0.0)
        self.assertAlmostEqual(system.cm_position[2], 0.0)

    def test_fold(self):
        system = copy.copy(self.ref)
        pos = copy.copy(system.particle[0].position)
        system.particle[0].position += system.cell.side
        system.particle[0].fold(system.cell)
        self.assertAlmostEqual(pos[0], system.particle[0].position[0])
        self.assertAlmostEqual(pos[1], system.particle[0].position[1])
        self.assertAlmostEqual(pos[2], system.particle[0].position[2])

    def test_overlaps(self):
        from atooms.system.particle import overlaps
        system = copy.copy(self.ref)
        for p in system.particle:
            p.radius = 1e-10
        pos = copy.copy(system.particle[1].position)
        system.particle[0].position = pos
        ov, ipart = overlaps(system.particle, system.cell)
        self.assertTrue(ov)
        self.assertEqual(ipart, [(0, 1)])

    def test_dump(self):
        self.assertEqual(self.ref.dump('spe')[-1],
                         self.ref.dump('particle.species')[-1])
        self.assertAlmostEqual(self.ref.dump('pos')[-1][-1],
                               self.ref.dump('particle.position')[-1][-1])
        self.assertAlmostEqual(self.ref.dump('vel')[-1][-1],
                               self.ref.dump('particle.velocity')[-1][-1])

    def test_species(self):
        system = copy.copy(self.ref)
        npart = len(system.particle)
        for p in system.particle[0: 10]:
            p.species = 'B'
        for p in system.particle[10: 30]:
            p.species = 'C'
        from atooms.system.particle import composition, distinct_species
        self.assertEqual(distinct_species(system.particle), ['A', 'B', 'C'])
        self.assertEqual(system.distinct_species(), ['A', 'B', 'C'])
        self.assertEqual(composition(system.particle)['A'], npart - 30)
        self.assertEqual(composition(system.particle)['B'], 10)
        self.assertEqual(composition(system.particle)['C'], 20)

    def test_packing(self):
        import math
        system = copy.copy(self.ref)
        self.assertAlmostEqual(system.packing_fraction * 6 / math.pi, system.density)

    def test_gyration(self):
        from atooms.system.particle import gyration_radius
        system = copy.copy(self.ref)

        # Ignore cell
        rg1 = gyration_radius(system.particle, method='N1')
        rg2 = gyration_radius(system.particle, method='N2')
        self.assertAlmostEqual(rg1, rg2)

        # With PBC all estimates are different but bounds must be ok
        rg1 = gyration_radius(system.particle, system.cell, method='min')
        rg2 = gyration_radius(system.particle, system.cell, method='N1')
        rg3 = gyration_radius(system.particle, system.cell, method='N2')
        self.assertLessEqual(rg1, rg2)
        self.assertLessEqual(rg3, rg2)

        # Equilateral triangle
        system.particle = [Particle(), Particle(), Particle()]
        system.particle[0].position = numpy.array([0.0, 0.0, 0.0])
        system.particle[1].position = numpy.array([1.0, 0.0, 0.0])
        system.particle[2].position = numpy.array([0.5, 0.5*3**0.5, 0])
        # Put the triangle across the cell
        system.particle[0].position -= 1.01*system.cell.side/2
        system.particle[1].position -= 1.01*system.cell.side/2
        system.particle[2].position -= 1.01*system.cell.side/2
        system.particle[0].fold(system.cell)
        system.particle[1].fold(system.cell)
        system.particle[2].fold(system.cell)
        rg1 = gyration_radius(system.particle, system.cell, method='min')
        rg2 = gyration_radius(system.particle, system.cell, method='N1')
        rg3 = gyration_radius(system.particle, system.cell, method='N2')
        self.assertAlmostEqual(rg1, 0.57735026919)
        self.assertAlmostEqual(rg2, 0.57735026919)
        self.assertAlmostEqual(rg3, 0.57735026919)

    def test_interaction(self):
        from atooms.interaction import Interaction
        system = copy.copy(self.ref)
        self.assertAlmostEqual(system.potential_energy(), 0.0)
        system.interaction = Interaction([])
        system.interaction.compute('energy', system.particle, system.cell)
        system.interaction.compute('forces', system.particle, system.cell)
        system.interaction.compute('stress', system.particle, system.cell)
        self.assertAlmostEqual(system.potential_energy(), 0.0)
        self.assertAlmostEqual(system.potential_energy(normed=True), 0.0)
        self.assertAlmostEqual(system.total_energy(), system.kinetic_energy())

    def test_overlap(self):
        from atooms.system.particle import self_overlap, collective_overlap
        sys1 = copy.deepcopy(self.ref)
        sys2 = copy.deepcopy(self.ref)
        sys1.particle = sys1.particle[:int(len(sys1.particle) / 2)]
        sys2.particle = sys2.particle[int(len(sys2.particle) / 2):]
        self.assertEqual(0, self_overlap(sys1.particle, sys2.particle, 0.001))
        self.assertEqual(0, collective_overlap(sys1.particle, sys2.particle, 0.001, sys1.cell.side))
        sys1.particle = sys1.particle
        sys2.particle = sys1.particle
        self.assertEqual(1, self_overlap(sys1.particle, sys2.particle, 0.001))
        self.assertEqual(1, collective_overlap(sys1.particle, sys2.particle, 0.001, sys1.cell.side))

    def test_overlap_random(self):
        # This test may fail from time to time
        from atooms.system.particle import collective_overlap
        N = 1000
        L = 5.0
        sys = [System(), System()]
        sys[0].cell = Cell([L, L, L])
        sys[1].cell = Cell([L, L, L])
        sys[0].particle = []
        sys[1].particle = []
        for _ in range(N):
            pos = [(random.random() - 0.5) * L,
                   (random.random() - 0.5) * L,
                   (random.random() - 0.5) * L]
            sys[0].particle.append(Particle(position=pos))
        for _ in range(N):
            pos = [(random.random() - 0.5) * L,
                   (random.random() - 0.5) * L,
                   (random.random() - 0.5) * L]
            sys[1].particle.append(Particle(position=pos))
        a = 0.3
        q_rand = ((a**3 * 4./3*3.1415) * N / sys[0].cell.volume)
        self.assertTrue(abs(q_rand - collective_overlap(sys[0].particle, sys[1].particle, a, sys[0].cell.side)) < 0.5)


if __name__ == '__main__':
    unittest.main()
