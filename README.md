Atooms
======

Atooms is a framework for classical, particle-based simulations written in python. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. This is achieved by offloading the critical parts of the calculation to backends written in C, CUDA or Fortran.

Quick start
-----------

Accessing the coordinates of a particle from an xyz trajectory file goes like this
```python
from atooms.trajectory import Trajectory

for system in Trajectory('input.xyz'):
    print system.particle[0].position
```

Let us change the density of the final configuration to unity
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    system.density = 1.0
```

We create a new trajectory file containing this "rescaled" configuration
```python
from atooms.trajectory import TrajectoryXYZ

with TrajectoryXYZ('output.xyz', 'w') as trajectory:
    trajectory.write(system, step=0)
```

We are ready to start a simulation from the rescaled configuration. We must specify which backend will carry out the simulation and pass it to a Simulation instance
```python
from atooms.backends.dryrun import DryRunBackend
from atooms.simulation import Simulation

backend = DryRunBackend(system)
sim = Simulation(backend)
sim.run(steps=10000)
print sim.system.temperature, sim.system.density
```
The DryRunBackend is just a dummy backend and won't do any actual simulation. Check out the available backends in the `backend` package or write your own!


Trajectory conversion
---------------------
Atooms ships with a command line tool to convert between various trajectory formats. For instance, the following will convert a trajectory file produced by rumd into a simpler xyz format

```bash
$ trj.py -i rumd -o xyz --stdout input.xyz.gz > output.xyz
```

`trj.py` has various options to control and format the output file. For instance, it is possible to include the particles' velocities in the output files by changing the output fields

```bash
$ trj.py -i rumd -o xyz --fmt-fields 'id,pos,vel' input.xyz.gz 
```

It is easy to add your own formats by subclassing any of the
Trajectory classes in the trajectory package. Just copy your custom
trajectory format as a python module under `atooms/plugins` and it
will become automatically available. 

Suppose your trajectory class is called `TrajectoryTest` and lies in
module `atooms/plugins/test.py`. You can convert a simple xyz file to
your new trajectory format like this

```bash
$ trj.py -i xyz -o test output.xyz > output.test
```

Features
--------
- High-level access to simulation objects
- Multiple trajectory formats
- Generic simulation interface with callback logic
- Efficient simulation backends, e.g. RUMD


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

Additional packages 
-------------------
atooms is composable: it is easy to add new functionalities, and just those you actually need.
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

Add the package root folder to $PYTHONPATH. Now you can import your package as

```python
import atooms.your_package
```

Authors
-------
Daniele Coslovich: http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/
