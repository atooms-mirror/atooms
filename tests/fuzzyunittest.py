from unittest import *

class FuzzyTestCase(TestCase):

    def assertAlmostEqual(x, y, reltol=1e-3):
        """Return True if x and y are equal to provided relative tolerance"""
        delta = abs(x-y)
        return self.assertLess(delta, reltol * 0.5*(x+y))

