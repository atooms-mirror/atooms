#!/usr/bin/env python

import os
import unittest
from atooms.core import utils


class Test(unittest.TestCase):

    def setUp(self):
        self.dirbase = '/tmp/test_core'
        self.dirout = self.dirbase + '/dir'
        utils.mkdir(self.dirout)
        with open(self.dirout + '/1', 'w'):
            pass
        with open(self.dirout + '/2', 'w'):
            pass
        utils.mkdir(self.dirout + '/3')

    def test_rm(self):
        utils.rmf([self.dirout + '/1', self.dirout + '/2'])
        self.assertFalse(os.path.exists(self.dirout + '/1'))
        self.assertFalse(os.path.exists(self.dirout + '/2'))
        self.assertTrue(os.path.exists(self.dirout + '/3'))

    def test_rm_glob(self):
        utils.rmf(self.dirout + '/*')
        self.assertFalse(os.path.exists(self.dirout + '/1'))
        self.assertFalse(os.path.exists(self.dirout + '/2'))
        self.assertTrue(os.path.exists(self.dirout + '/3'))

    def test_rmd(self):
        utils.rmd(self.dirbase)


if __name__ == '__main__':
    unittest.main(verbosity=0)
