Atooms
======

[![pypi](https://img.shields.io/pypi/v/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![version](https://img.shields.io/pypi/pyversions/atooms.svg)](https://pypi.python.org/pypi/atooms/)
[![license](https://img.shields.io/pypi/l/atooms.svg)](https://en.wikipedia.org/wiki/GNU_General_Public_License)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1183301.svg)](https://doi.org/10.5281/zenodo.1183301)

**atooms** is a Python framework for simulations of interacting particles. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. To achieve this, atooms relies on backends written in C, CUDA or Fortran.

Quick start
-----------

The goal of atooms is to provide a coherent interface to the basic objects of particle simulations, such as [molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) or [Monte Carlo](https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_physics) simulations. 
The simulation data are usually stored in trajectory files, which atooms makes it easy to analyze, manipulate and convert.

In this simple example, we read a trajectory file in [xyz format](https://en.wikipedia.org/wiki/XYZ_format) composed by multiple frames. Each frame holds the state of the system at a given instant of time during the simulation. Accessing the coordinates of the particles in a trajectory file goes like this:
```python
from atooms.trajectory import Trajectory

with Trajectory('input.xyz') as trajectory:
    for system in trajectory:
        print('The position of particle 0 is', system.particle[0].position)
```
Note that trajectories support iteration and indexing, just like lists.

Here we pick the last frame of the trajectory, change the density of the system to unity and write this new configuration to a trajectory format suitable for the [RUMD](http://rumd.org) simulation package:
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    system.density = 1.0
    print('The new density is', len(system.particle) / system.cell.volume)

from atooms.trajectory import TrajectoryRUMD
with TrajectoryRUMD('rescaled.xyz.gz', 'w') as trajectory:
    trajectory.write(system, step=0)
```

We can now run 1000 molecular dynamics steps using the Lennard-Jones potential:
```python
from atooms.backends.rumd import RUMD
from atooms.simulation import Simulation

backend = RUMD('rescaled.xyz.gz', forcefield_file='lj_rumd.ff', 
               output_path='/tmp/outdir', integrator='nve')
sim = Simulation(backend)
sim.run(1000)
print('Final temperature and density', sim.system.temperature, sim.system.density)
```
The forcefield file `lj_rumd.ff` (available in `data/`) defines the interaction potential.

Documentation
-------------
See the [tutorial](https://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/atooms/) for a step-by-step introduction to atooms objects and the [public API documentation](https://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/docs/atooms/) for full details. 

Installation
------------
From the python package index
```
pip install atooms
```

From the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms.git
cd atooms
make install
```

Simulation backends
-------------------
atooms has a generic simulation interface that abstracts out most of the common parts of particle-based simulations. The actual simulation code is wrapped by a simulation backend that exposes a minimal but consistent interface. This enables one to develop more complex simulation frameworks (e.g., [parallel tempering](https://gitlab.info-ufr.univ-montp2.fr/atooms/parallel_tempering)) that are essentially decoupled from the underlying simulation code.

Trajectory conversion
---------------------
atooms provides a command line tool to convert between various trajectory formats. The following command will convert a trajectory file produced by [RUMD](http://rumd.org) into a simpler xyz format

```bash
$ trj.py convert -i rumd -o xyz trajectory.xyz.gz output.xyz
```
If you don't specify the output path, the trajectory is written to standard output. This is useful for quick inspection of complex trajectory formats or for piping into sed / awk.

`trj.py` provides means to fine tune the format of the output file. Type `trj.py --help` to get a list of options and supported trajectory formats.

Custom trajectory formats 
-------------------------
It is easy to add new trajectory formats by subclassing existing trajectory classes. Just create a package called
`atooms_plugins` and add your trajectory modules there. They will be automatically
available to all client codes that use atooms.

Suppose you wrote a custom trajectory class `TrajectoryABC` in
`atooms_plugins/test.py` (the last path is relative to the current
directory). You can now convert an existing xyz trajectory to your custom
format:

```bash
$ trj.py convert output.xyz output.abc
```

Remember to add an empty `__init__.py` file at the root of `atooms_plugins`. 
Actually, the `atooms_plugins` package can be put anywhere in your `PYTHONPATH`.

Additional packages 
-------------------
Atooms is composable: it is easy to add new functionalities, and just those you actually need.
Additional packages are available from the [atooms main repository](https://gitlab.info-ufr.univ-montp2.fr/atooms).
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
Daniele Coslovich: http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/
