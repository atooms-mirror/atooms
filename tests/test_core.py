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
        utils.cp(self.dirout + '/1', self.dirout + '/2')
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

    def test_wget(self):
        utils.wget('https://framagit.org/atooms/atooms/-/raw/master/README.md', self.dirbase)
        utils.download('https://framagit.org/atooms/atooms/-/raw/master/README.md', self.dirbase)

    def test_timer(self):
        t = utils.Timer()
        t.start()
        t.stop()
        t.wall_time
        t.cpu_time

    def test_tipify(self):
        from atooms.core.utils import tipify
        self.assertTrue(type(tipify("2.0")) is float)
        self.assertTrue(type(tipify("2")) is int)
        self.assertTrue(type(tipify("t2")) is str)


if __name__ == '__main__':
    unittest.main()
