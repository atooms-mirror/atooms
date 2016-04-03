#!/usr/bin/env python

import sys
from atooms import trajectory

trj_map = {
    'auto': trajectory.Trajectory,
    'xyz' : trajectory.TrajectoryXYZ,
    'rumd': trajectory.TrajectoryRUMD,
    'hoomd': trajectory.TrajectoryHOOMD,
    'lammps': trajectory.TrajectoryLAMMPS,
    'hdf5': trajectory.TrajectoryHDF5,
    }

def change_density_keep_N(s, args):
    rho_old = s.density
    x = (rho_old / args['rho'])**(1./3)
    s.cell.side *= x
    for p in s.particle:
        p.position *= x
    print '# changed density: %s -> %s' % (rho_old, s.density)
    return s

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--density', dest='density', action='store', default=None, type=float, help='new density')
    parser.add_argument(      '--step',    dest='step', action='store', default=None, type=int, help='')
    parser.add_argument('-i', '--fmt-inp', dest='inp', type=str, default='auto', help='input format ')
    parser.add_argument('-o', '--fmt-out', dest='out', type=str, default='', help='output format for conversion')
    parser.add_argument(nargs='+', dest='files',type=str, help='input files')
    args = parser.parse_args()

    for f in args.files:
        if args.density is None:
            continue
        with trj_map[args.inp](f) as t:
            fout = trajectory.convert(t, t.__class__, tag='-rho%s' % args.density, exclude=('vx','vy','vz'),
                                      callback=change_density_keep_N, args={'rho':args.density})
    
