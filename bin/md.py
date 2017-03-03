#!/usr/bin/env python

"""Simple molecular dynamics simulation driver using RUMD backend."""

import os
from atooms.simulation import Simulation, Scheduler, WriterThermo, WriterConfig
from atooms.backends.rumd_backend import RumdBackend
from atooms.utils import setup_logging

def main(params):
    # TODO: dump params to a file in output_dir
    if params.verbose:
        setup_logging(level=10)
    with open(params.forcefield) as f:
        exec(f.read())
    potential = locals()['_potential']
    if params.T is not None:
        params.integrator = 'nvt'
    s = RumdBackend(params.input_file, potential,
                    integrator=params.integrator,
                    temperature=params.T, dt=params.dt)
    sa = Simulation(s, output_path=os.path.join(params.output_dir, 'trajectory'),
                     checkpoint_interval=params.config_interval,
                     steps=params.steps,
                     restart=params.restart)
    sa.add(WriterThermo(), Scheduler(params.thermo_interval))
    sa.add(WriterConfig(), Scheduler(params.config_interval))
    sa.run()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-T',   dest='T', type=float, help='temperature')
    parser.add_argument('--dt', dest='dt', type=float, default=0.002, help='time step')
    parser.add_argument('-i',   dest='input_file', help='input_file')
    parser.add_argument('-I',   dest='integrator', default='nve', help='integrator')
    parser.add_argument('--ff', dest='forcefield', help='force field')
    parser.add_argument('-n',   dest='steps', type=int, default=0, help='number of steps')
    parser.add_argument('-t','--thermo-interval', dest='thermo_interval', type=int, default=0, help='energy interval')
    parser.add_argument('-c','--config-interval', dest='config_interval', type=int, default=0, help='config interval')
    parser.add_argument('-r',   dest='restart', action='store_true', help='restart')
    parser.add_argument('-v',   dest='verbose', action='store_true', help='verbose output')
    parser.add_argument(dest='output_dir',type=str, help='output directory')
    params = parser.parse_args()

    main(params)
