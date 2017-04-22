#!/usr/bin/env python

import unittest
from atooms.trajectory import utils

class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_block_period(self):
        tpl = [0, 1, 2, 4, 8, 16, 32]
        more = tpl[:-1]
        more += [tpl[-1]+i for i in tpl[:-1]]
        self.assertEqual(utils.get_period(more), 7)
        self.assertEqual(utils.check_block_period(more, utils.get_period(more)), None)

if __name__ == '__main__':
    unittest.main(verbosity=0)


