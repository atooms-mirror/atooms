#!/usr/bin/env python

"""Simple molecular dynamics simulation driver using RUMD backend."""

import sys
import os
from atooms.core import __version__, __commit__
from atooms.simulation import Simulation, Scheduler, write_thermo, write_config
from atooms.simulation.backends import RumdBackend
from atooms.utils import setup_logging, report_parameters, report_command, mkdir, cp

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
    cp(params.forcefield, output_base + '.ff')
    cp(params.forcefield, output_base + '.chk.ff')
    report_parameters(params.__dict__, output_base + '.params', '%s+%s' % (__version__, __commit__))
    report_command(sys.argv[0], params.__dict__, ['output_dir'], output_base + '.cmd')

    s = RumdBackend(params.input_file, params.forcefield,
                    integrator=params.integrator,
                    temperature=params.T, dt=params.dt,
                    fixcm_interval=params.fixcm_interval)

    sa = Simulation(s, output_path=output_base,
                     checkpoint_interval=params.checkpoint_interval,
                     steps=params.steps,
                     restart=params.restart)

    if params.backend_output:
        s._suppress_all_output = False
        s._initialize_output = True
        s.rumd_simulation.sample.SetOutputDirectory(output_base)
        s.rumd_simulation.SetOutputScheduling("energies", "linear", interval=params.thermo_interval)
        s.rumd_simulation.SetOutputScheduling("trajectory", params.config_sampling, interval=params.config_interval)
        s.rumd_simulation.SetOutputMetaData("trajectory", precision=6, virials=False)
        if params.config_interval > 0 and params.config_sampling == "logarithmic":
            s.rumd_simulation.SetBlockSize(params.config_interval)
        else:
            s.rumd_simulation.SetBlockSize(params.steps)
        # Trim target steps to be a multiple of config_interval
        # params.steps = params.steps / params.config_interval * params.config_interval
    else:
        sa.add(write_thermo, Scheduler(params.thermo_interval))
        sa.add(write_config, Scheduler(params.config_interval))
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
    parser.add_argument(     '--config-sampling', dest='config_sampling', default="logarithmic", help='')
    parser.add_argument('-t','--thermo-interval', dest='thermo_interval', type=int, default=0, help='energy interval')
    parser.add_argument('-c','--config-interval', dest='config_interval', type=int, default=0, help='config interval')
    parser.add_argument('-C','--checkpoint-interval', dest='checkpoint_interval', type=int, default=0, help='checkpoint interval')
    parser.add_argument(     '--fixcm-interval', dest='fixcm_interval', type=int, default=100, help='fix cm interval')
    parser.add_argument('-r',   dest='restart', action='store_true', help='restart')
    parser.add_argument('-v',   dest='verbose', action='store_true', help='verbose output')
    parser.add_argument('-d',   dest='debug', action='store_true', help='debug output')
    parser.add_argument('-b',   dest='backend_output', action='store_true', help='use backend output')
    parser.add_argument(dest='output_dir',type=str, help='output directory')
    params = parser.parse_args()

    main(params)
