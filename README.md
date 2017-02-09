Atooms
======

Atooms is a python framework for particle-based simulations. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. This is achieved by offloading the critical parts of the calculation to backends written in C, CUDA or Fortran.

Quick start
-----------

This simple example will show how to manipulate particles from a trajectory file. Accessing the coordinates of the first particle in your system goes like this
```python
from atooms.trajectory import Trajectory

for system in Trajectory('input.xyz'):
    print 'The position of particle 0 is', system.particle[0].position
```
In this example, the trajectory format (xyz) will be simply guessed from the file extension.

Let us change the density of the final configuration to unity
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    system.density = 1.0
```
Note that trajectories composed by multiple samples support iteration and slicing, just like lists.

Now we create a new trajectory file containing this "rescaled" configuration
```python
from atooms.trajectory import TrajectoryXYZ

with TrajectoryXYZ('output.xyz', 'w') as trajectory:
    trajectory.write(system, step=0)
```
Here, we explicitly used the `TrajectoryXYZ` class instead of the factory class `Trajectory`.

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
Atooms ships with a command line tool to convert between various trajectory formats. For instance, the following will convert a trajectory file produced by [RUMD](http://rumd.org) into a simpler xyz format

```bash
$ trj.py -i rumd -o xyz input.xyz.gz output.xyz
```
Dropping the output file, will dump the new trajectory to standard output.

`trj.py` has various options to control the format of the output file. For instance, it is possible to include the particles' velocities in the output file by changing the output fields

```bash
$ trj.py -i rumd -o xyz --fmt-fields 'id,pos,vel' input.xyz.gz output.xyz
```
To get a list of the available trajectory formats

```bash
$ trj.py --fmt-available
```

Plugins
-------
It is easy to add a new trajectory format by subclassing any of the
existing Trajectory classes. Just create a package called
`atooms_plugins` and add your trajectory modules there. They will be automatically
available to all client codes that use atooms.

Suppose you added a custom trajectory class `TrajectoryABC` to
`atooms_plugins/test.py` (the last path is relative to the current
directory). You can convert an existing xyz trajectory to your custom
format like this

```bash
$ trj.py output.xyz output.abc
```

Remember to add an empty `__init__.py` file at the root of `atooms_plugins`. 
Actually, the `atooms_plugins` package can be put anywhere in your `PYTHONPATH`.

Additional features
-------------------
- High-level access to simulation objects
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

Add the package root folder to $PYTHONPATH. You can now import your package as

```python
import atooms.your_package
```

Authors
-------
Daniele Coslovich: http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/
