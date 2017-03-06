Atooms
======

Atooms is a python framework for particle-based simulations. It makes it easy to develop simulation and analysis tools using an expressive language, without sacrificing efficiency. This is achieved by offloading the critical parts of the calculation to backends written in C, CUDA or Fortran.

Quick start
-----------

This simple example will show how to manipulate particles from a trajectory file and run a simulation using one of the molecular dynamics backends. Accessing the coordinates of the particles in your system goes like this:
```python
from atooms.trajectory import Trajectory

for system in Trajectory('input.xyz'):
    print 'The position of particle 0 is', system.particle[0].position
```

Let us change the density of the final configuration to unity:
```python
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    system.density = 1.0
```
Note how trajectories with multiple frames support iteration and slicing, just like lists.

We store this new configuration in a trajectory file suitable for RUMD:
```python
from atooms.trajectory import TrajectoryRUMD

with TrajectoryRUMD('rescaled.xyz.gz', 'w') as trajectory:
    trajectory.write(system, step=0)
```

We are ready to start an NVE simulation from the rescaled configuration:

```python
from atooms.backends.rumd import RumdBackend
from atooms.simulation import Simulation

backend = RumdBackend('rescaled.xyz.gz', output_path='/tmp/output_dir', integrator='nve')
sim = Simulation(backend)
sim.run(steps=10000)
print sim.system.temperature, sim.system.density
```

Trajectory conversion
---------------------
Atooms provides a command line tool to convert between various trajectory formats. The following command will convert a trajectory file produced by [RUMD](http://rumd.org) into a simpler xyz format

```bash
$ trj.py -i rumd -o xyz input.xyz.gz output.xyz
```
If you don't specify the output path, the trajectory is written to standard output. This useful for quick inspection of complex trajectory formats or for piping into sed / awk.

`trj.py` has various options to control the format of the output file. For instance, it is possible to include the particles' velocities in the output file by changing the output fields:

```bash
$ trj.py -i rumd -o xyz --fmt-fields 'id,pos,vel' input.xyz.gz output.xyz
```
Type `trj.py --help` to get a list of options and supported trajectory formats.


Custom trajectory formats 
-------------------------
It is easy to add a new trajectory format by subclassing any of the
existing trajectory classes. Just create a package called
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
