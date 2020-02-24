#!/usr/bin/env python

"""Test trajectory decorators."""

import os
import unittest
import numpy

from atooms.trajectory import TrajectoryXYZ


def _equal(system1, system2):
    for p1, p2 in zip(system1.particle, system2.particle):
        if p1.species != p2.species:
            print(p1.species, p2.species)
            return False
    return True


class TestDecorators(unittest.TestCase):

    Trajectory = TrajectoryXYZ

    def setUp(self):
        self.finp = '/tmp/test_decorators.xyz'
        with open(self.finp, 'w') as fh:
            fh.write("""\
2
step:1 cell:6.0,6.0,6.0
A 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
step:2 cell:6.0,6.0,6.0
A 1.1 -1.1 0.0
A -2.9 -2.9 0.0
2
step:3 cell:6.0,6.0,6.0
A 1.2 -1.2 0.0
A -2.9 2.9 0.0
2
step:4 cell:6.0,6.0,6.0
A 1.2 -1.2 0.0
A -2.9 2.9 0.0
""")

    def test_unfolded_jump(self):
        from atooms.trajectory.decorators import Unfolded
        with self.Trajectory(self.finp) as th:
            with Unfolded(th) as tu:
                self.assertEqual(list(tu[2].particle[1].position), [3.1, -3.1, 0.0])  # unfolded
                self.assertEqual(list(th[2].particle[1].position), [-2.9, 2.9, 0.0])  # original

    def test_sliced(self):
        from atooms.trajectory.decorators import Sliced
        with Sliced(self.Trajectory(self.finp), slice(0, 2, 1)) as th:
            self.assertEqual(th.steps, [1, 2])
            self.assertEqual(list(th[-1].particle[0].position), [1.1, -1.1, 0.0])
            self.assertEqual(list(th[-1].particle[1].position), [-2.9, -2.9, 0.0])

    def test_callback_args(self):
        def cbk_1(system, scale):
            for p in system.particle:
                p.position *= scale
            return system

        def cbk_2(system, offset):
            for p in system.particle:
                p.position -= offset
            return system

        def cbk_3(system, memory):
            """Callback that modifies a mutable."""
            memory.append(1)
            return system

        memo = []
        with self.Trajectory(self.finp) as th:
            th.register_callback(cbk_1, 2.0)
            th.register_callback(cbk_2, offset=1.0)
            th.register_callback(cbk_3, memo)
            self.assertEqual(th.steps, [1, 2, 3, 4])
            self.assertEqual(list(th[-1].particle[0].position),
                             [1.2 * 2 - 1.0, -1.2 * 2 - 1.0, 0.0 * 2 - 1.0])
            self.assertEqual(memo, [1])

    def test_change_species(self):
        from copy import deepcopy
        from atooms.system import System, Particle
        system_A = System([Particle(species='A'), Particle(species='B')])
        system_C = System([Particle(species='0'), Particle(species='1')])
        system_F = System([Particle(species='1'), Particle(species='2')])
        from atooms.trajectory.decorators import change_species
        # DO nothing here
        self.assertTrue(_equal(system_A, change_species(deepcopy(system_A), 'A')))
        self.assertTrue(_equal(system_C, change_species(deepcopy(system_C), 'C')))
        self.assertTrue(_equal(system_F, change_species(deepcopy(system_F), 'F')))
        # Change
        self.assertTrue(_equal(system_C, change_species(deepcopy(system_A), 'C')))
        self.assertTrue(_equal(system_F, change_species(deepcopy(system_A), 'F')))
        self.assertTrue(_equal(system_A, change_species(deepcopy(system_C), 'A')))
        self.assertTrue(_equal(system_F, change_species(deepcopy(system_C), 'F')))
        self.assertTrue(_equal(system_A, change_species(deepcopy(system_F), 'A')))
        self.assertTrue(_equal(system_C, change_species(deepcopy(system_F), 'C')))

    def tearDown(self):
        os.remove(self.finp)


if __name__ == '__main__':
    unittest.main()
