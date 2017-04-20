#!/usr/bin/env python

"""Simple molecular dynamics simulation driver using RUMD backend."""

import sys
import os
from atooms.core import __version__, __commit__
from atooms.simulation import Simulation, Scheduler, WriterThermo, WriterConfig
from atooms.simulation.backends import RumdBackend
from atooms.utils import setup_logging, report_parameters, report_command, mkdir

def main(params):
    if params.verbose:
        setup_logging(level=20)
    if params.debug:
        setup_logging(level=10)
    if params.T is not None:
        params.integrator = 'nvt'
    if os.path.exists(params.input_file + '.ff'):
        params.forcefield = params.input_file + '.ff'
    output_base = os.path.join(params.output_dir, 'trajectory')
    mkdir(output_base)
    report_parameters(params.__dict__, output_base + '.params', '%s+%s' % (__version__, __commit__))
    report_command(sys.argv[0], params.__dict__, ['output_dir'], output_base + '.cmd')
    s = RumdBackend(params.input_file, params.forcefield,
                    integrator=params.integrator,
                    temperature=params.T, dt=params.dt,
                    fixcm_interval=params.fixcm_interval)
    sa = Simulation(s, output_path=output_base,
                     checkpoint_interval=params.config_interval,
                     steps=params.steps,
                     restart=params.restart)

    if params.backend_output:
        s._suppress_all_output = False
        s._initialize_output = True
        s.rumd_simulation.sample.SetOutputDirectory(output_base)
        s.rumd_simulation.SetOutputScheduling("energies", "linear", interval=params.thermo_interval)
        s.rumd_simulation.SetOutputScheduling("trajectory", "logarithmic")
        if params.config_interval > 0:
            s.rumd_simulation.SetBlockSize(params.config_interval)
        else:
            s.rumd_simulation.SetBlockSize(params.steps)
        # Trim target steps to be a multiple of config_interval
        # params.steps = params.steps / params.config_interval * params.config_interval
    else:
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
    parser.add_argument(     '--fixcm-interval', dest='fixcm_interval', type=int, default=1000, help='fix cm interval')
    parser.add_argument('-r',   dest='restart', action='store_true', help='restart')
    parser.add_argument('-v',   dest='verbose', action='store_true', help='verbose output')
    parser.add_argument('-d',   dest='debug', action='store_true', help='debug output')
    parser.add_argument('-b',   dest='backend_output', action='store_true', help='use backend output')
    parser.add_argument(dest='output_dir',type=str, help='output directory')
    params = parser.parse_args()

    main(params)
