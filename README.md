Atooms
======

[![pypi](https://img.shields.io/pypi/v/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![version](https://img.shields.io/pypi/pyversions/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![license](https://img.shields.io/pypi/l/atooms.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183301.svg)](https://doi.org/10.5281/zenodo.1183301)
[![pipeline](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)
[![coverage report](https://framagit.org/atooms/atooms/badges/master/coverage.svg)](https://framagit.org/atooms/atooms/-/commits/master)

**atooms** is a Python framework for simulations of interacting particles. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. To achieve this, atooms relies on backends written in C, CUDA or Fortran.

Quick start
-----------

The goal of atooms is to provide a coherent interface to the basic objects of [molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) or [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_physics) simulations.
The simulation data are stored in trajectory files, which are easy to analyze, manipulate and convert with atooms.

Here we read a multi-frame trajectory in [xyz format](https://en.wikipedia.org/wiki/XYZ_format). Each frame holds the state of the system at a given instant of time during the simulation. Accessing the coordinates of a particle goes like this:
```python
from atooms.trajectory import Trajectory

with Trajectory('input.xyz') as trajectory:
    for system in trajectory:
        print(system.particle[0].position)
```
Note that trajectories support iteration and indexing, just like lists.

Here we pick the last frame of the trajectory, change the density of the system to unity and write this new configuration to a trajectory format suitable for the [RUMD](http://rumd.org) simulation package:
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    system.density = 1.0
    print('New density:', len(system.particle) / system.cell.volume)

from atooms.trajectory import TrajectoryRUMD
with TrajectoryRUMD('rescaled.xyz.gz', 'w') as trajectory:
    trajectory.write(system)
```

atooms has a generic simulation interface that abstracts out most of the common parts of particle-based simulations. The actual simulation code is wrapped by a simulation backend that exposes a minimal but consistent interface. This enables one to develop more complex simulation frameworks (e.g., [parallel tempering](https://framagit.org/atooms/parallel_tempering)) that are essentially decoupled from the underlying simulation code.

Here we can now a short molecular dynamics simulation with the RUMD backend, using a Lennard-Jones potential:
```python
import rumd
from atooms.backends.rumd import RUMD
from atooms.simulation import Simulation

potential = rumd.Pot_LJ_12_6(cutoff_method=rumd.ShiftedPotential)
potential.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
backend = RUMD('rescaled.xyz.gz', [potential], integrator='nve'
sim = Simulation(backend)
sim.run(1000)
print('Final temperature and density:', sim.system.temperature, sim.system.density)
```

Documentation
-------------
Check out the [tutorial](https://atooms.frama.io/atooms/tutorial) for more examples and the [public API](https://atooms.frama.io/api/atooms) for full details.

Installation
------------
From the python package index
```
pip install atooms
```

From the code repository
```
git clone https://framagit.org/atooms/atooms.git
cd atooms
make install
```

Contributing
------------
You are welcome to contribute to this project! Please have a look at [these guidelines](https://framagit.org/atooms/atooms/-/blob/atooms-3.0.0/CONTRIBUTING.md).

Additional packages 
-------------------
Atooms is composable: it is easy to add new functionalities, and just those you actually need.
Additional packages are available from the [atooms main repository](https://framagit.org/atooms).
These packages will be installed in the atooms namespace to prevent name clashing.

If you want to add your own package to the atooms namespace, structure it this way
```bash
atooms/your_package
atooms/your_package/__init__.py
```

where ```__init__.py``` contains

```python
from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)
```

Add the package root folder to $PYTHONPATH. You can now import your package as

```python
import atooms.your_package
```

Authors
-------
Daniele Coslovich: https://www.units.it/daniele.coslovich/
