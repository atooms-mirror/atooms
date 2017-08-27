#!/usr/bin/env python
"""
Correlate particles properties from trajectory files.

Example:
--------

  tpaste.py file1.xyz:radius file2.xyz.voronoi.xyz:volume
"""

import sys
from atooms import trajectory as trj

f1, attr1 = sys.argv[1].split(':')
f2, attr2 = sys.argv[2].split(':')
t1 = trj.Trajectory(f1)
t2 = trj.Trajectory(f2)

for step, s1, s2 in trj.utils.paste(t1, t2):
    try:
        for i in range(len(s1.particle)):
            print getattr(s1.particle[i], attr1), getattr(s2.particle[i], attr2)
    except:
        print getattr(s1, attr1), getattr(s2, attr2)
