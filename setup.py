#!/usr/bin/env python

import os
import subprocess
from distutils.core import setup

setup(name='atooms',
      #version=git_version,
      description='Python tools for atomistic simulations',
      author='Daniele Coslovich',
      author_email='daniele.coslovich@univ-montp2.fr',
      #http://www.coulomb.univ-montp2.fr/perso/daniele.coslovich/
      url='',
      packages=[],
      # 'atooms', 'atooms.adapters', 'atooms.system',
      # 'atooms.interaction', 'atooms.potential', 
      # 'atooms.postprocessing'],
      scripts=['bin/trj.py', 'bin/md.py', 'bin/modes.py', 'bin/vis.py'],
      #scripts=glob.glob(os.path.join('bin', '*.py'))
     )
