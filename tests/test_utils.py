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
B 1.0 -1.0 0.0
1
step:3 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
""")
        with TrajectoryXYZ(finp) as th:
            self.assertTrue(utils.is_cell_variable(th, tests=-1))
            self.assertFalse(utils.is_cell_variable(th, tests=1))
            self.assertTrue(utils.is_cell_variable(th, tests=2))

            # Test utility function
            from atooms.trajectory.utils import dump
            pos = dump(th, what='pos')
            from atooms.trajectory.utils import is_semigrandcanonical, is_grandcanonical
            self.assertFalse(is_grandcanonical(th))
            self.assertTrue(is_semigrandcanonical(th, tests=3))

    def test_formats(self):
        from atooms.trajectory.utils import formats
        self.assertTrue('xyz' in formats())

    def test_file_index(self):
        from atooms.trajectory import utils
        finp = '/tmp/test_utils.xyz'
        with open(finp, 'w') as fh:
            fh.write("""\
2
step:0 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
A 1.0 -1.0 0.0
3
step:1 cell:2.0,1.0,1.0
A 1.0 -1.0 0.0
A 1.0 -1.0 0.0
A 1.0 -1.0 0.0
3
step:2 cell:1.0,1.0,1.0
A 1.0 -1.0 0.0
A 1.0 -1.0 0.0
A 1.0 -1.0 0.0
""")
        with open(finp) as fh:
            utils.file_index(fh)

    def test_dumps(self):
        from atooms.core import utils
        utils.report_command('cmd', {}, '', None)
        utils.report_parameters({'x': 1.0}, None, '1.0.0')

    def test_timings(self):
        from atooms.core import utils
        f = lambda x : x
        utils.clockit(f)
        f(1)
        t = utils.Timer()
        t.start()
        t.stop()
        t.cpu_time
        t.wall_time

    def tearDown(self):
        from atooms.core.utils import rmf
        rmf('/tmp/test_utils*')


if __name__ == '__main__':
    unittest.main()
