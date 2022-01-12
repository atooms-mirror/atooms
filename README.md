Atooms
======

[![pypi](https://img.shields.io/pypi/v/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![version](https://img.shields.io/pypi/pyversions/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![license](https://img.shields.io/pypi/l/atooms.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fatooms/HEAD?labpath=docs%2F)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183301.svg)](https://doi.org/10.5281/zenodo.1183301)
[![pipeline](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)
[![coverage report](https://framagit.org/atooms/atooms/badges/master/coverage.svg)](https://framagit.org/atooms/atooms/-/commits/master)

**atooms** is a Python framework for simulations of interacting particles. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. To achieve this, atooms relies on backends written in C, CUDA or Fortran.

Quick start
-----------

The goal of the base atooms package is to provide a coherent interface to the basic objects of [molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) or [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_physics) simulations.

As an example, we set up a mixture of two types of particles in an elongated box with rough walls on the x-axis boundaries formed by a third component. We set the number density to 1 (in reduced units).
```python
from atooms.system import System

system = System(N=64)
system.replicate(times=4, axis=0)
system.density = 1.0
system.composition = {'A': len(system.particle) - 40, 'B': 40}
for p in system.particle:
    if abs(p.position[0]) > 0.7 * system.cell.side[0] / 2:
        # Freeze particles on the cell borders along the x-axis
        p.mass = 1e100
        p.species = 'C'
    else:
        # Randomly displace those in the bulk
        import numpy.random
        p.position += (numpy.random.random() - 0.5) * 0.3
        p.fold(system.cell)
system.show('ovito')
```

![](https://framagit.org/atooms/atooms/-/raw/master/snapshot.png)

Simulation data are stored in trajectory files, which are easy to manipulate and convert with atooms. Here, we write our system in a single-frame trajectory file using the [xyz format](https://en.wikipedia.org/wiki/XYZ_format).

```python
from atooms.trajectory import TrajectoryXYZ

with TrajectoryXYZ('input.xyz', 'w') as th:
    th.variables = ['species', 'position', 'mass']
    th.write(system)
```

The trajectory file can be used to start a simulation using one of the supported [simulation backends](https://atooms.frama.io/atooms/tutorial/simulations.html) or your own code.

Documentation
-------------
Check out the [tutorial](https://atooms.frama.io/atooms/tutorial) for more examples and the [public API](https://atooms.frama.io/api/atooms) for full details.

Org-mode and jupyter notebooks are available under `docs/`. You can run the tutorial interactively on [Binder]( https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fatooms/HEAD?labpath=docs%2).

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

Component packages 
------------------
Atooms is modular: it is easy to add new functionalities, and just those you actually need.
Component packages are available from the [atooms main repository](https://framagit.org/atooms).

Here a few things one can do with them:

- compute [static and dynamic correlation functions](https://framagit.org/atooms/postprocessing)
- use the [models repository](https://framagit.org/atooms/models) to build an interaction model for simple liquids and glasses
- explore and analyze the [potential energy landscape](https://framagit.org/atooms/landscape)
- run efficient [multi-GPU parallel tempering]((https://framagit.org/atooms/parallel_tempering) simulations

Component packages are installed in the atooms namespace to prevent name clashing. If you want to add your own package to the atooms namespace, structure it this way
```bash
atooms/your_package
atooms/your_package/__init__.py
```

where ```__init__.py``` contains

```python
from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)
```

Install `your_package` and you are ready to go
```python
import atooms.your_package
```

Authors
-------
Daniele Coslovich: https://www.units.it/daniele.coslovich/
