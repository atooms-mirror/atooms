Atooms
======

[![pypi](https://img.shields.io/pypi/v/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![version](https://img.shields.io/pypi/pyversions/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![license](https://img.shields.io/pypi/l/atooms.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fatooms/HEAD?labpath=docs%2F)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183301.svg)](https://doi.org/10.5281/zenodo.1183301)
[![pipeline](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)
[![coverage report](https://framagit.org/atooms/atooms/badges/master/coverage.svg)](https://framagit.org/atooms/atooms/-/commits/master)

**atooms** is a Python framework for simulations of interacting particles. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency.

This is the base atooms package, which provides some core classes and functionalities. For instance, it allows you to
- create [starting configurations](https://atooms.frama.io/atooms/tutorial/basics.html) for your simulations
- convert between [trajectory file](https://atooms.frama.io/atooms/tutorial/trajectories.html) formats - even in-house ones
- run molecular dynamics simulations using one of the [supported backends](https://atooms.frama.io/atooms/tutorial/simulations.html)

Several [component packages](Component packages) are built on top of the core library. Here a few things one can do with them:
- compute [static and dynamic correlation functions](https://framagit.org/atooms/postprocessing)
- use the [models repository](https://framagit.org/atooms/models) to build an interaction model for simple liquids and glasses
- explore and analyze the [potential energy landscape](https://framagit.org/atooms/landscape)
- run efficient [multi-GPU parallel tempering]((https://framagit.org/atooms/parallel_tempering) simulations

Quick start
-----------

The goal of the base atooms package is to provide a coherent interface to the basic objects of [molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) or [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_physics) simulations.

As an example, we set up a mixture of two types of particles, A and B, in an elongated box. The number density is set to unity.
```python
import numpy
from atooms.system import System

system = System(N=64)
system.replicate(times=4, axis=0)
system.composition = {'A': 128, 'B': 128}
system.density = 1.0
```

Particles in the central part of the cell get a random displacement.
```python
for p in system.particle:
    # Randomly displace particles in the bulk
    if abs(p.position[0]) < system.cell.side[0] / 4:
        p.position += 0.5 * (numpy.random.random() - 0.5)
        p.fold(system.cell)
system.show('ovito')
```

![](https://framagit.org/atooms/atooms/-/raw/master/snapshot.png)

Simulation data are stored in trajectory files, which are easy to manipulate and convert with atooms. Here, we write the system species and positions in a single-frame trajectory file using the [xyz format](https://en.wikipedia.org/wiki/XYZ_format).
```python
from atooms.trajectory import TrajectoryXYZ

with TrajectoryXYZ('input.xyz', 'w') as th:
    th.variables = ['species', 'position']  # actually, this is the default
    th.write(system)
```

This trajectory file can be used to start a simulation using one the available [simulation backends](https://atooms.frama.io/atooms/tutorial/simulations.html) or your own code.

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

They are installed in the atooms namespace to prevent name clashing. If you want to add your own package to the atooms namespace, structure it this way
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
