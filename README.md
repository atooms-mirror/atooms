Atooms
======

[![pypi](https://img.shields.io/pypi/v/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![version](https://img.shields.io/pypi/pyversions/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![license](https://img.shields.io/pypi/l/atooms.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fatooms/HEAD?labpath=docs%2F)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183301.svg)](https://doi.org/10.5281/zenodo.1183301)
[![pipeline](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)](https://framagit.org/atooms/atooms/badges/master/pipeline.svg)
[![coverage report](https://framagit.org/atooms/atooms/badges/master/coverage.svg)](https://framagit.org/atooms/atooms/-/commits/master)

**atooms** is a high-level Python framework for simulations of interacting particles, such as molecular dynamics or Monte Carlo.

This is the core package: it provides a consistent interface to the basic objects of particle-based simulations. [Feature packages](Component packages) are built on top of it and implement complex simulation methods and analysis tools.

Quick start
-----------

Here is a small example for setting up a mixture of two types of particles, A and B, in a periodic elongated cell. The number density is set to unity.
```python
from atooms.system import System

system = System(N=64)
system.replicate(times=4, axis=0)
system.composition = {'A': 128, 'B': 128}
system.density = 1.0
```

Particles in the central part of the cell get a random displacement and are folded back into the simulation cell
```python
import numpy

for p in system.particle:
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

The trajectory file can now be used to start a simulation using one the available [simulation backends](https://atooms.frama.io/atooms/tutorial/simulations.html) or your own code.

Features
--------

- Focus on a simple and expressive interface
- API refined over the years towards consistency
- Modular and extensible design via namespace packages
- Semantic versioning - for what is worth
- Easy to interface: in-house codes and custom formats are first-class citizens
- Support for efficient simulation backends, with a focus on GPU codes

Documentation
-------------
Check out the [tutorial](https://atooms.frama.io/atooms/tutorial) for more examples and the [public API](https://atooms.frama.io/atooms/api/atooms) for more details.

Org-mode and jupyter notebooks are available under [docs](https://framagit.org/atooms/atooms/-/blob/master/docs/). You can run them interactively on [Binder](https://mybinder.org/v2/git/https%3A%2F%2Fframagit.org%2Fatooms%2Fatooms/HEAD?labpath=docs%2).

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
You are welcome to contribute to this project! Please have a look at [these guidelines](https://framagit.org/atooms/atooms/-/blob/master/CONTRIBUTING.md).

Feature packages 
------------------
Atooms is modular: it is easy to add new functionalities, and just those you actually need.

Feature packages are available from the [atooms main repository](https://framagit.org/atooms). They are installed in the `atooms` namespace to prevent name clashing. If you want to add your own feature package to the atooms namespace, structure it this way
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

Thanks go to Francesco Turci and Geert Kapteijns for their contributions.