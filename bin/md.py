#!/usr/bin/env python

"""Molecular dynamics simulation driver"""

import sys
import os
from atooms.adapters.rumd_backend import Simulation, single
from atooms.simulation import log
from atooms.simulation.parallel_tempering import ParallelTempering

def main(params):
    # TODO: dump params to a file in output_dir
    if params.verbose:
        log.setLevel(10)
    execfile(params.forcefield)
    potential = locals()['_potential']
    s = single(params.input_file,
               potential,
               params.T,
               params.dt)
    sa = Simulation(s, fixcm_interval=params.fixcm_interval,
                     output_path=os.path.join(params.output_dir, 'configs'), 
                     thermo_interval=params.thermo_interval,
                     config_interval=params.config_interval,
                     checkpoint_interval=params.config_interval,
                     steps=params.steps,
                     restart=params.restart)
    sa.run()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='verbose', action='store_true', help='verbose output')
    parser.add_argument('-T', dest='T', type=float, help='temperature')
    parser.add_argument('--dt', dest='dt', type=float, default=0.004, help='time step')
    parser.add_argument('-t','--thermo-interval', dest='thermo_interval', type=int, default=50000, help='config interval (in units of swap periods)')
    parser.add_argument('-c','--config-interval', dest='config_interval', type=int, default=50000, help='energy interval (in units of swap periods)')
    parser.add_argument('--fixcm-interval', dest='fixcm_interval', type=int, default=100, help='interval in steps between CM fixing)')
    parser.add_argument('-n','--steps', dest='steps', type=int, default=None, help='number of steps')
    parser.add_argument('-i', dest='input_file', help='input_file')
    parser.add_argument('--ff', dest='forcefield', help='force field')
    parser.add_argument('-r', dest='restart', action='store_true', help='restart')
    parser.add_argument(dest='output_dir',type=str, help='output directory')
    params = parser.parse_args()

    main(params)
