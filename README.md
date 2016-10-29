Atooms
======

This is atooms - a framework for atomistic simulations. It enables development of simulation and analysis tools using a natural, expressive language but without sacrificing efficiency. This is achieved by splitting code into high-level python layers and backends written in C, CUDA or Fortran.

Quick start
-----------

Accessing the coordinates of a particle from a multi-configuration xyz trajectory file goes like this
```python
from atooms.trajectory import Trajectory

for system in Trajectory('input.xyz'):
    print system.particle[0].position
```

Let us compress the final configuration of the trajectory to unit density
```python
rho = 1.0
with Trajectory('input.xyz') as trajectory:
    system = trajectory[-1]
    factor = (system.density / rho)**(1./3)
    for particle in system.particle:
        particle.position *= factor
    system.cell.side *= factor
```
Actually, the System class has a ```rescale()``` method to do just that. 

We create a new trajectory file with just the rescaled configuration
```python
from atooms.trajectory import TrajectoryXYZ

with TrajectoryXYZ('output.xyz', 'w') as trajectory:
    trajectory.write(system, step=0)
```

To start a simulation from the rescaled configuration, we pick up a simulation backend and pass it to a Simulation object
```python
from atooms.backends.dryrun import DryRunBackend
from atooms.simulation import Simulation

backend = DryRunBackend(system)
sim = Simulation(backend)
sim.run(steps=10000)
print sim.system.temperature, sim.system.density
```
Of course, we just pretended to do 10000 steps! The DryRunBackend won't do any actual simulation, nor write anything to disk. Check out the available backends or write your own!


Features
--------
- High-level access to simulation objects
- Multiple trajectory formats
- Generic simulation interface with callback logic
- Efficient simulation backends, e.g. RUMD

Additional packages 
-------------------
atooms is composable: it makes it easy to add new functionalities and just those you actually need.
Some additional packages are available here https://gitlab.info-ufr.univ-montp2.fr/atooms.
Once installed, they will be accessible in the atooms namespace.

If you want to add your package to the atooms namespace, structure it this way
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
Daniele Coslovich <daniele.coslovich@umontpellier.fr>
