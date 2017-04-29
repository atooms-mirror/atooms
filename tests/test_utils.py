#!/usr/bin/env python

import unittest
from atooms.trajectory import utils

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
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), None)

    def test_block_size_3(self):
        more = [0, 1, 10, 12, 20, 30]
        self.assertEqual(utils.get_block_size(more), len(more))

    def test_block_size_3(self):
        more =  [0, 1, 2, 4, 8, 16]
        more += [32, 33, 34, 36, 40, 48]
        self.assertEqual(utils.get_block_size(more), 6)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), None)

    def test_block_size_4(self):
        more =  [0, 1, 2, 4, 8, 16]
        more += [32 + i for i in more]
        self.assertEqual(utils.get_block_size(more), 6)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), None)

    def test_block_size_5(self):
        more =  [0, 1, 2, 4, 8, 16, 24, 32]
        more += [40 + i for i in more]
        self.assertEqual(utils.get_block_size(more), 8)
        self.assertEqual(utils.check_block_size(more, utils.get_block_size(more)), None)

if __name__ == '__main__':
    unittest.main(verbosity=0)


