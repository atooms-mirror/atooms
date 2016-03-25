#!/usr/bin/env python

import sys
from vapory import *


_particle = []
_bonds = []
_radius = {"A": 0.2, "B": 0.15}
_color = {"A": (1,1,1), "B": (0.5,0.5,1.0)}

def box(cell, radius=0.1):
    L = cell.side[0]/2
    return [Cylinder((-L,-L,-L), (L,-L,-L), radius, 'color', (0,0,0)),
            Cylinder((-L,-L,-L), (-L,L,-L), radius),
            Cylinder((-L,-L,-L), (-L,-L,L), radius),
            Cylinder((-L,-L,L), (L,-L,L), radius),
            Cylinder((-L,-L,L), (-L,L,L), radius),
            Cylinder((-L,L,L), (-L,L,-L), radius),
            Cylinder((-L,L,L), (L,L,L), radius),
            Cylinder((-L,L,-L), (L,L,-L), radius),
            Cylinder((L,L,-L), (L,L,L), radius),
            Cylinder((L,L,-L), (L,-L,-L), radius),
            Cylinder((L,-L,-L), (L,-L,L), radius),
            Cylinder((L,-L,L), (L,L,L), radius),
    ]

def show(particle, radius=None, color=None, rscale=1.0, rcut=0.0, append=False):
    if radius is None:
        radius = [_radius[p.name.strip()] for p in particle]
    if color is None:
        color = [_color[p.name.strip()] for p in particle]
    for p, r, c in zip(particle, radius, color):
        _particle.append(Sphere(p.position, r*rscale if r>rcut else 0, Texture(Pigment('color', c),
                                                                               Finish('phong', 0.0,
                                                                                      'reflection', 0.0))))
    return _particle

def main(f, args):
    from atooms.trajectory import Trajectory
    t = Trajectory(f)
    s = t[0]
    L = s.cell.side[0]

    particles = show(s.particle, radius=[p.radius for p in s.particle], rscale=1.)
    lines = box(s.cell)
    camera = Camera('location', [1.5*L, 0.7*L, 0.8*L], 'look_at', [3,-2,0] )
    light = [LightSource([11*L, 3*L, 0],  'color', [.9,.9,.9], 'shadowless'),
             LightSource([11*L, -3*L, 0],  'color', [.9,.9,.9], 'shadowless')
    ]
    scene = Scene(camera, objects= particles + light + lines + [Background((1,1,1))])
    if args.output is None:
        scene.render(f + ".step%d.png" % 0, width=800, height=800)
    else:
        scene.render(args.output, width=800, height=800)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-v', dest='verbose', action='store_true', help='verbose output')
parser.add_argument('-o', dest='output', action='store', default=None, help='output file')
parser.add_argument('--rc', dest='rcut', action='store', default=None, type=float, help='output file')
parser.add_argument(nargs='+', dest='files',type=str, help='input files')
args = parser.parse_args()

for f in args.files:
    main(f, args)
    
