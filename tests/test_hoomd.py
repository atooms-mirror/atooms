#!/usr/bin/env python

import unittest

from atooms import trajectory

class Test(unittest.TestCase):

    def setUp(self):
        self.finp = '/tmp/test.xyz'
        fh = open(self.finp, 'w')
        fh.write("""\
2
1
A 1.0 -1.0 0.0
A 2.9 -2.9 0.0
2
2
A 1.1 -1.1 0.0
A -2.9 -2.9 0.0
2
3
A 1.2 -1.2 0.0
A -2.9 2.9 0.0
2
4
A 1.3 -1.3 0.0
A -2.8 2.8 0.0
6.0 6.0 6.0
""")
        fh.close()

    def test(self):
        fname = '/tmp/test.tgz'
        txyz = trajectory.TrajectoryXYZ(self.finp)

        t = trajectory.TrajectoryHOOMDXML(fname, 'w:gz')
        t.write_sample(txyz.read_sample(0), 0, 0)
        t.close()
        t = trajectory.TrajectoryHOOMDXML(fname, 'r')
        
        self.assertEqual(t.read_sample(0).particle[0].position[1], 
                             txyz.read_sample(0).particle[0].position[1])

        self.assertEqual(t.read_sample(0).cell.side[1], 
                             txyz.read_sample(0).cell.side[1])

        t.close()

        txyz.close()

if __name__ == '__main__':
    unittest.main()


