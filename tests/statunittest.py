"""Unit test enhancements

unittest in python 2.7 has a convenient argument delta in
assertAlmostEqual. However, this is missing in python 2.6. To avoid
requiring unittest2, which backports new features introduced in 2.7,
we provide our custom subclass of TestCase.
"""

import unittest

class StatTestCase(unittest.TestCase):
    
    def assertAlmostEqual(self, first, second, places=7, msg=None, delta=None):
        if delta is None:
            unittest.TestCase.assertAlmostEqual(self, first, second, places, msg)
        else:
            return abs(first-second) <= delta
