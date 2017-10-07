#!/usr/bin/env python

import os
import shutil
import unittest
import numpy

from atooms.core.utils import mkdir, rmd
from atooms.trajectory.folder import TrajectoryFolder


class TestXYZ(unittest.TestCase):

    def setUp(self):
        self.dirname = '/tmp/test_folder'
        rmd(self.dirname)
        mkdir(self.dirname)
        for i in range(10, 13):
            fname = os.path.join(self.dirname, 'step_%d' % i)
            with open(fname, 'w') as fh:
                fh.write("""\
2
step:%d
A 1.0 -1.0 0.0
A 2.9 -2.9 0.0
""" % i)

    def test_xyz_indexed(self):
        r_ref = [[1., -1., 0.], [2.9, -2.9, 0.]]
        t = TrajectoryFolder(self.dirname, step_pattern='step_(\d*)')
        self.assertEqual(t.steps, [10, 11, 12])

    def tearDown(self):
        if os.path.exists(self.dirname):
            shutil.rmtree(self.dirname)


if __name__ == '__main__':
    unittest.main()
