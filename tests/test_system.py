#!/usr/bin/env python

import random
import copy
import numpy
import unittest
from atooms.system import System
from atooms.system.cell import Cell
from atooms.system.particle import Particle

class Test(unittest.TestCase):

    def setUp(self):
        N = 100
        L = 10.0
        self.ref = System()
        self.ref.cell = Cell([L, L, L])
        self.ref.particle = []
        while len(self.ref.particle) <= N:
            pos = [(random.random()-0.5) * L,
                   (random.random()-0.5) * L,
                   (random.random()-0.5) * L]
            self.ref.particle.append(Particle(position=pos))

    def test_density(self):
        system = copy.copy(self.ref)
        density_old = system.density
        system.density = density_old*1.1
        self.assertAlmostEqual(system.density, density_old*1.1)

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
            p.radius = 0.1
        pos = copy.copy(system.particle[1].position)
        ov, ipart = overlaps(system.particle, system.cell)
        self.assertFalse(ov)
        system.particle[0].position = pos
        ov, ipart = overlaps(system.particle, system.cell)
        self.assertTrue(ov)
        self.assertEqual(ipart, [(0, 1)])

    def test_dump(self):
        self.assertAlmostEqual(self.ref.dump('ids')[-1],
                               self.ref.dump('particle.id')[-1])
        self.assertAlmostEqual(self.ref.dump('pos')[-1][-1],
                               self.ref.dump('particle.position')[-1][-1])
        self.assertAlmostEqual(self.ref.dump('vel')[-1][-1],
                               self.ref.dump('particle.velocity')[-1][-1])

if __name__ == '__main__':
    unittest.main()


