#!/usr/bin/env python

import unittest
from atooms.trajectory import utils
from atooms.trajectory import TrajectoryXYZ


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_block_size_1(self):
        more = [0, 1, 2, 3, 4]
        self.assertEqual(utils.get_block_size(more), 1)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), None)

    def test_block_size_2(self):
        more = [0, 1, 10, 11, 20, 21]
        self.assertEqual(utils.get_block_size(more), 2)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), more)

    def test_block_size_3(self):
        more = [0, 1, 10, 12, 20, 30]
        self.assertEqual(utils.get_block_size(more), len(more))

    def test_block_size_3(self):
        more = [0, 1, 2, 4, 8, 16]
        more += [32, 33, 34, 36, 40, 48]
        self.assertEqual(utils.get_block_size(more), 6)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), more)

    def test_block_size_4(self):
        more = [0, 1, 2, 4, 8, 16]
        more += [32 + i for i in more]
        self.assertEqual(utils.get_block_size(more), 6)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), more)

    def test_block_size_5(self):
        more = [0, 1, 2, 4, 8, 16, 24, 32]
        more += [40 + i for i in more]
        self.assertEqual(utils.get_block_size(more), 8)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), more)

    def test_cell_variable(self):
        finp = '/tmp/test_utils.xyz'
#         with open(finp, 'w') as fh:
#             fh.write("""\
# 1
# step:0 cell:1.0,1.0,1.0
# A 1.0 -1.0 0.0
# 1
# step:1 cell:1.0,1.0,1.0
# A 1.0 -1.0 0.0
# 1
# step:2 cell:1.0,1.0,1.0
# A 1.0 -1.0 0.0
# 1
# step:3 cell:1.0,1.0,1.0
# A 1.0 -1.0 0.0
# """)
#         t = Trajectory(finp)
#         self.assertFalse(utils.is_cell_variable(t, tests=-1))

        with open(finp, 'w') as fh:
            fh.write("""\
1
step:0 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:1 cell:2.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:2 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
1
step:3 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
""")
        with TrajectoryXYZ(finp) as th:
            self.assertTrue(utils.is_cell_variable(th, tests=-1))
            self.assertFalse(utils.is_cell_variable(th, tests=1))
            self.assertTrue(utils.is_cell_variable(th, tests=2))

    def tearDown(self):
        from atooms.core.utils import rmf
        rmf('/tmp/test_utils*')


if __name__ == '__main__':
    unittest.main()
