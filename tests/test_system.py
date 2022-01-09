#!/usr/bin/env python

import sys
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
        random.seed(1)
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

    def test_ndim(self):
        system = System()
        self.assertEqual(system.number_of_dimensions, 0)
        system.cell = Cell([2.0, 2.0])
        self.assertEqual(system.number_of_dimensions, 2)

    def test_density(self):
        system = copy.copy(self.ref)
        density_old = system.density
        system.density = density_old * 1.1
        self.assertAlmostEqual(system.density, density_old * 1.1)

        # When there are no particles, density is zero
        system = System()
        self.assertAlmostEqual(system.density, 0.0)

    def test_density_cluster(self):
        system = copy.deepcopy(self.ref)
        system.cell = None
        self.assertAlmostEqual(system.density, 0.1, 1)

    def test_temperature(self):
        system = copy.copy(self.ref)
        system.set_temperature(1.0)
        self.assertAlmostEqual(system.temperature, 1.0)
        system.scale_velocities(1.0)
        self.assertAlmostEqual(system.temperature, 1.0)

        # Pathological case
        system.particle = system.particle[0: 1]
        self.assertAlmostEqual(system.temperature, 2.0 / system.number_of_dimensions * system.kinetic_energy())

        # Empty system
        system = System()
        self.assertAlmostEqual(system.temperature, 0.0)

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

    def test_pbc_center(self):
        system = copy.copy(self.ref)
        # Move the center of the cell so that positions are within 0 and L
        system.cell.center = system.cell.side / 2
        for p in system.particle:
            p.position += system.cell.side / 2
        # Check that distances are the same
        for i in range(len(system.particle)):
            self.assertAlmostEqual(sum(self.ref.particle[0].distance(self.ref.particle[i], self.ref.cell)**2),
                                   sum(system.particle[0].distance(system.particle[i], system.cell)**2))
        # Move one particle out of the box and fold it back
        pos = copy.copy(system.particle[0].position)
        system.particle[0].position += system.cell.side
        system.particle[0].fold(system.cell)
        self.assertAlmostEqual(pos[0], system.particle[0].position[0])
        self.assertAlmostEqual(pos[1], system.particle[0].position[1])
        self.assertAlmostEqual(pos[2], system.particle[0].position[2])

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
        self.assertEqual(system.distinct_species, ['A', 'B', 'C'])
        self.assertEqual(composition(system.particle)['A'], npart - 30)
        self.assertEqual(composition(system.particle)['B'], 10)
        self.assertEqual(composition(system.particle)['C'], 20)

        # Unhashable numpy scalars
        for p in system.particle:
            p.species = numpy.array(1)
        self.assertEqual(system.distinct_species, [1])

    def test_species_layout(self):
        system = copy.copy(self.ref)
        for p in system.particle[0: 10]:
            p.species = 'B'
        for p in system.particle[10: 30]:
            p.species = 'C'
        self.assertTrue(system.species_layout == 'A')
        system.species_layout = 'C'
        self.assertTrue(system.species_layout == 'C')
        system.species_layout = 'F'
        self.assertTrue(system.species_layout == 'F')
        system.species_layout = 'A'
        self.assertTrue(system.species_layout == 'A')

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
        from atooms.system.interaction import InteractionBase
        system = copy.copy(self.ref)
        self.assertAlmostEqual(system.potential_energy(), 0.0)
        system.interaction = InteractionBase()
        system.interaction.compute('energy', system.dump('position'))
        self.assertAlmostEqual(system.potential_energy(), 0.0)
        self.assertAlmostEqual(system.potential_energy(normed=True), 0.0)
        self.assertAlmostEqual(system.total_energy(), system.kinetic_energy())

    def test_interaction_add(self):
        from atooms.system.interaction import InteractionBase
        x, y = InteractionBase(), InteractionBase()
        x.energy, y.energy = 1., 1.
        z = sum([x, y])
        self.assertAlmostEqual(z.energy, 2.)

        # This must fail
        ok = False
        y.energy = None
        try:
            z = x + y
        except:
            ok = True
        self.assertTrue(ok)

    def test_interaction_terms(self):
        from atooms.system.interaction import InteractionBase, Interaction
        x, y = InteractionBase(), InteractionBase()
        z = Interaction(x, y)
        w = Interaction()
        w.add(x)
        w.add(y)
        w.compute('forces', position=numpy.array([[1.]]))
        self.assertAlmostEqual(w.energy, 0.0)
        self.assertTrue(w.forces is not None)

    def test_interaction_terms_via_system(self):
        from atooms.system.interaction import InteractionBase, Interaction
        system = copy.copy(self.ref)
        system.interaction = Interaction(InteractionBase())
        system.compute_interaction('forces')
        self.assertAlmostEqual(system.potential_energy(), 0.0)

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

    def test_view(self):
        import numpy
        from atooms.system import Particle, System

        p = [Particle(), Particle()]
        s = System(p)
        pos = s.dump("pos", order='F', view=True)

        # Modify the dumped array in place preserves the view
        pos[:, 0] += 1.0
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))
        # Modify the position array in place preserves the view
        p[0].position *= 2
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))
        # Modify the position array in place preserves the view
        p[0].position[:] = p[0].position[:] + 4
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))
        pos[:, 0] = pos[:, 0] + 1.0
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))
        # Reassining the position will of course destroy the view
        p[0].position = p[0].position * 2
        self.assertFalse((p[0].position == pos[:, 0]).all())
        self.assertFalse(numpy.may_share_memory(p[0].position, pos[:, 0]))

    def test_view_clear(self):
        import numpy
        from atooms.system import Particle, System

        p = [Particle(), Particle()]
        s = System(p)

        # We check that particle positions are views on dump array
        pos = s.dump("pos", order='F', view=True)
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))

        # We should get the same dump array
        pos = s.dump("pos", order='F', view=True)
        self.assertTrue((p[0].position == pos[:, 0]).all())
        self.assertTrue(numpy.may_share_memory(p[0].position, pos[:, 0]))

        # We clear the dump array
        pos1 = s.dump("pos", order='F', view=True, clear=True)
        pos1 += 1
        self.assertFalse((p[0].position == pos[:, 0]).all())
        self.assertFalse(numpy.may_share_memory(p[0].position, pos[:, 0]))

    def test_dump_cbk(self):
        """
        Make sure that applying callbacks or changing the particle array
        size creates a new dump.
        """
        import numpy
        from atooms.system import Particle, System

        for view in [False, True]:
            p = [Particle(), Particle()]
            s = System(p)
            pos1 = s.dump("pos", view=view)

            def cbk(system):
                s = copy.copy(system)
                s.particle = [system.particle[0]]
                return s
            s = cbk(s)
            pos2 = s.dump('pos', view=view)
            self.assertEqual(pos1.shape, (2, 3))
            self.assertEqual(pos2.shape, (1, 3))
            # Grandcanonical
            s.particle.append(Particle())
            self.assertEqual(s.dump('pos', view=view).shape, (2, 3))
            # Reassign particle
            # Expected failure with view = True
            if not view:
                s.particle[0] = Particle(position=[1.0, 1.0, 1.0])
                self.assertEqual(s.dump('pos', view=view)[0][0], 1.0)

    def test_dump_species(self):
        """
        Make sure that changing species in the dump is reflected in the
        particle species and viceversa.
        """
        import numpy
        from atooms.system import Particle, System

        view = True
        p = [Particle(), Particle()]
        s = System(p)
        spe = s.dump("particle.species", view=True)
        spe[0] = 'B'
        self.assertEqual(spe[0], s.particle[0].species)
        # With this syntax, the numpy scalar preserves the view!
        # We should upgrade species to property and hide this inside
        s.particle[0].species[()] = 'C'
        self.assertEqual(spe[0], s.particle[0].species)

    def test_dumps(self):
        """Check that dump order does not matter"""
        import numpy
        from atooms.system import Particle, System

        p = [Particle(), Particle()]
        s = System(p)
        # View, dump, view: we should preserve the view and get a copy of the dump
        pos1 = s.dump("pos", view=True)
        pos1[0, 0] = 0.0
        pos2 = s.dump('pos', view=False)
        pos2[0, 0] = 1.0
        pos3 = s.dump("pos", view=True)
        self.assertNotAlmostEqual(pos1[0, 0], pos2[0, 0])
        self.assertAlmostEqual(pos1[0, 0], pos3[0, 0])
        self.assertFalse(numpy.may_share_memory(pos1, pos2))
        self.assertTrue(numpy.may_share_memory(pos1, pos3))

    def test_dump_flatten(self):
        """Check that flattening a dump does not change successive views"""
        import numpy
        from atooms.system import Particle, System

        p = [Particle(), Particle()]
        s = System(p)
        self.assertEqual(len(s.dump('pos', view=True).shape), 2)
        self.assertEqual(len(s.dump('pos', view=False, flat=True).shape), 1)
        self.assertEqual(len(s.dump('pos', view=True).shape), 2)
        self.assertEqual(len(s.dump('pos', view=False, flat=True).shape), 1)

    def test_dump_fail(self):
        """Check that dump fails on unknown attributes"""
        for what in ['whatever', 'particle.whatever', 'cell.whatever']:
            try:
                self.ref.dump(what)
            except AttributeError:
                pass

    def test_decimate(self):
        from atooms.system import Particle, System
        from atooms.system.particle import composition, decimate
        p = [Particle(species='A')]*20 + [Particle(species='B')]*10
        pnew = decimate(p, 12)
        x = composition(pnew)
        self.assertEqual(x['A'], 8)
        self.assertEqual(x['B'], 4)

    @unittest.skipIf(sys.version_info.major == 2, 'skip show() tests with python 2')
    def test_show(self):
        N = 3
        L = 5.0
        system = System()
        system.cell = Cell([L, L, L])
        system.particle = []
        for _ in range(N):
            pos = (numpy.random.random(len(system.cell.side)) - 0.5) * system.cell.side
            p = Particle(position=pos)
            system.particle.append(p)

        try:
            system.show(backend='')
        except ValueError:
            pass
        try:
            system.show(backend='matplotlib')
            system.show(backend='ovito')
            system.show(backend='3dmol')
        except ImportError:
            self.skipTest('missing backend')

    def test_rotate(self):
        """Rotate particles so that the principal axis is along the y axis"""
        from atooms.system import Particle, Cell
        from atooms.system.particle import rotate
        p1, p2, p3 = Particle(position=[0.0, 0.0, 0.0]), Particle(position=[1.0, 0.0, 0.0]), Particle(position=[2.0, 0.0, 0.0])
        particle = [p1, p2, p3]
        cell = Cell(side=[10.0, 10.0, 10.0])
        rotated = rotate(particle, cell)
        self.assertAlmostEqual(rotated[0].position[1], 0, 6)
        self.assertAlmostEqual(rotated[1].position[1], 1, 6)
        self.assertAlmostEqual(rotated[2].position[1], 2, 6)

    def test_init(self):
        """Test __init__() syntax"""
        self.assertEqual(len(System(N=1).particle), 1)
        self.assertEqual(len(System(N=10).particle), 10)
        self.assertEqual(len(System(N=100).particle), 100)
        self.assertEqual(System(N=10).composition, {'A': 10})
        self.assertEqual(System(N={'B': 4}).composition, {'B': 4})
        self.assertEqual(System(N={'Si': 4, 'O': 2}).composition, {'Si': 4, 'O': 2})
        self.assertEqual(System(N={'A': 4, 'B': 10, 'C': 2}).composition, {'A': 4, 'B': 10, 'C': 2})
        self.assertAlmostEqual(System(N={'Si': 10, 'O': 20}).concentration, {'Si': 10/30., 'O': 20/30.})
            

if __name__ == '__main__':
    unittest.main()
