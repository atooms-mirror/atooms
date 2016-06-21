Atooms
======

This is atooms, a general-purpose, object-oriented package for molecular simulations.

Copyright (C) 2010-2016 Daniele Coslovich <daniele.coslovich@umontpellier.fr>

Plugins
-------

Put your plugin packages under local/ add the path to your
PYTHONPATH. In order for atooms to recognize your plugins under the
namespace package atooms.plugins, your plugin structure should be
the following

local/<yourplugin>/atooms/plugins
local/<yourplugin>/atooms/plugins/__init__.py

where __init__.py contains

```python
from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)
```