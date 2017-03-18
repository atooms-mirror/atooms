Atooms
======

Atooms is a python framework for particle-based simulations. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. This is achieved by offloading the critical parts of the calculation to backends written in C, CUDA or Fortran.

Quick start
-----------

The idea is to provide a convenient interface to the basic objects at the core of molecular dynamics or Monte Carlo simulations.
In this simple example, we read a trajectory file in [xyz format](https://en.wikipedia.org/wiki/XYZ_format). Accessing the coordinates of the particles in a trajectory file goes like this:
```python
from atooms.trajectory import Trajectory

with Trajectory('input.xyz') as trajectory:
    for system in trajectory:
        print 'The position of particle 0 is', system.particle[0].position
```
Note that trajectories support iteration and slicing, just like lists. 

We can modify the density of the final configuration to unity and store this new configuration in a different trajectory format. Here we use a format suitable for the [RUMD](http://rumd.org) backend:
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    factor = (system.density / rho)**(1./3)
    for particle in system.particle:
        particle.position *= factor
    system.cell.side *= factor

from atooms.trajectory import TrajectoryRUMD
with TrajectoryRUMD('trajectory.xyz.gz', 'w') as trajectory:
    trajectory.write(system, step=0)
```

Installation
------------
From the python package index [coming up soon!]
```
pip install atooms
```

Alternatively, from the code repository
```
git clone https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms.git
cd atooms
make install
```

Trajectory conversion
---------------------
Atooms provides a command line tool to convert between various trajectory formats. The following command will convert a trajectory file produced by [RUMD](http://rumd.org) into a simpler xyz format

```bash
$ trj.py -i rumd -o xyz trajectory.xyz.gz output.xyz
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
$ trj.py output.xyz output.abc
```

Remember to add an empty `__init__.py` file at the root of `atooms_plugins`. 
Actually, the `atooms_plugins` package can be put anywhere in your `PYTHONPATH`.

Simulation backends
-------------------

Atooms has a generic simulation interface that abstracts out most of the common parts of particle-based simulations. The actual simulation code is wrapped by a simulation "backend" that exposes a minimal but unified interface. This enables one to develop more complex simulation frameworks (e.g., [parallel tempering](https://gitlab.info-ufr.univ-montp2.fr/atooms/parallel_tempering)) that are essentially decoupled from the underlying simulation code.

This is a quick example how to run 10000 molecular dynamics steps using the [RUMD](http://rumd.org) backend:

```python
from atooms.backends.rumd import RumdBackend
from atooms.simulation import Simulation

backend = RumdBackend('rescaled.xyz.gz', forcefield_file='lj_rumd.py', 
                      output_path='/tmp/outdir', integrator='nve')
sim = Simulation(backend)
sim.run(10000)
print 'Final temperature and density', sim.system.temperature, sim.system.density
```
The forcefield file `lj_rumd.py` (available in `data/`) defines the interaction potential.

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
