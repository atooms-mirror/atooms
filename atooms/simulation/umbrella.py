class Umbrella(object):

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self, sim):
        print 'calling', self.func.__name__, self.args[0].__name__, self.kwargs
        return self.func(sim, *self.args, **self.kwargs)

def f(x):
    return x

#u = Umbrella(f, 1.0)

def quadratic_umbrella_len(sim, k, x_0):
    x = len(sim.trj)
    return 0.5 * k * (x - x_0)**2

def quadratic_umbrella(sim, obs, k, x_0):
    x = obs(sim)
    return 0.5 * k * (x - x_0)**2

def bias(sim, obs, s):
    x = obs(sim)
    return - s * x

# def obs(sim):
#     return len(sim.trj)

# def mobility(sim):
#     return mobility(sim.trj)

# import os
# import logging
# import numpy as np
# import atooms.core.progress
# from atooms.trajectory.decorators import Unfolded, filter_species
# from atooms.simulation import Simulation
# from atooms.system import Thermostat
# from atooms.backends.dryrun import DryRun
# from atooms.backends.lammps import LAMMPS
# from atooms.core.utils import setup_logging, mkdir
# from atooms.transition_path_sampling import core, TransitionPathSampling
# from atooms.trajectory import TrajectoryXYZ

# lmp = DryRun() #LAMMPS('', '')
# sim = Simulation(lmp, steps=1)
# tps = TransitionPathSampling(sim, 1.0)
# tps.umbrella = Umbrella(quadratic_umbrella, obs, k=1.0, x_0=1.0)
# print tps.umbrella(tps)
