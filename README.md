Atooms
======

This is atooms - a composable, high-level framework for molecular simulations.

Copyright (C) 2010-2016 Daniele Coslovich <daniele.coslovich@umontpellier.fr>

Examples
--------

Access particles' coordinates in a multi-configuration xyz trajectory file
```python
from atooms.trajectory import Trajectory

for system in Trajectory('input.xyz'):
    print system.particle[0].position
```

Features
--------
- High-level interface to simulation objects
- Handle and convert multiple trajectory formats 
- Flexible simulation interface with callback logic
- Efficient molecular dynamics backends, e.g. RUMD

Adding packages to atooms namespace
-----------------------------------

Structure your package this way

atooms/your_package
atooms/your_package/__init__.py

where __init__.py contains

```python
from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)
```

Add the package root folder to $PYTHONPATH. You can now import your package as

```python
import atooms.your_package
```

