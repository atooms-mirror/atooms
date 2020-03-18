"""
Generic Umbrella class for umbrella sampling and biased simulations
"""

# Default umbrella functions

def bias(sim, observable, s):
    x = observable(sim)
    return - s * x

def quadratic_umbrella(sim, observable, k, x_0):
    x = observable(sim)
    return 0.5 * k * (x - x_0)**2

def quadratic_umbrella_len(sim, k, x_0):
    x = len(sim.trj)
    return 0.5 * k * (x - x_0)**2


class Umbrella(object):

    """Callable class that represents a generic biasing potential"""

    def __init__(self, func, *args, **kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return '{} {}: {}'.format(self.func.__name__, self.args, self.kwargs)

    def __call__(self, sim):
        """The umbrella takes a Simulation instance as a first argument"""
        return self.func(sim, *self.args, **self.kwargs)
