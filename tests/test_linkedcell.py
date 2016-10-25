#!/usr/bin/env python

import unittest
import numpy
from atooms.system.linkedcell import LinkedCell

class TestLinkedCell(unittest.TestCase):

    def setUp(self):
        pass

    def test_(self):
        n = 10
        rc = numpy.array([2, 2, 2])
        box = numpy.array([10, 10, 10])
        pos = numpy.array([[1, 1, 1],
                           [1, 1, 5]])

        l = LinkedCell()
        l.adjust(n, box, rc)
        l.compute(pos, box)

if __name__ == '__main__':
    unittest.main(verbosity=0)


