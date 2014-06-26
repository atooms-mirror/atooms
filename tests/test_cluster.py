#!/usr/bin/env python

import unittest
from atooms.cluster import ClusterAnalysis

class TestCluster(unittest.TestCase):

    def setUp(self):
        pass

    def test_simple_voronoi(self):
        inp = [[(0,2,8),   (0,2,8,1)], 
               [(0,2,3,4), (0,2,3,4), (0,2,1)], 
               [(0,2,8),   (0,3,6)]]
        ref = [set([(0, 2, 8, 1), (0, 2, 8), (0, 3, 6)]), 
               set([(0, 2, 3, 4), (0, 2, 1)])]
        c = ClusterAnalysis(inp)
        self.assertEqual(c.clusters, ref)

#     def test_from_file_voronoi(self):
#         fh = open('/tmp/test_cluster.txt', 'w')
#         fh.write("""\
# 0 2 8
# 0 2 3 4
# 0 2 8

# 0 2 8 1
# 0 2 3 4 1
# 0 3 6 1

# 0 2 8 2
# 0 2 8
# 0 3 6
# """)
#         fh.close()
#         ref = [set([(0, 2, 8), (0, 2, 3, 4), (0, 2, 8, 2), (0, 3, 6)]), set([(0, 2, 8, 1), (0, 3, 6, 1), (0, 2, 3, 4, 1)])]
#         c = Cluster(fh.name)
#         self.assertEqual(c.cluster, ref)

    def test_from_file_index(self):
        fh = open('/tmp/test_cluster.txt', 'w')
        fh.write("""\
# n 2 t = 10
1 2 3 
2 8 10

# n 3 t = 20
1 2 4 
5 8 10
10 12
""")
        fh.close()
        #ref = [set([(0, 2, 8), (0, 2, 3, 4), (0, 2, 8, 2), (0, 3, 6)]), set([(0, 2, 8, 1), (0, 3, 6, 1), (0, 2, 3, 4, 1)])]
        c = ClusterAnalysis(fh.name)
        print c.cluster
        self.assertEqual(c.cluster, ref)

if __name__ == '__main__':
    unittest.main(verbosity=0)


