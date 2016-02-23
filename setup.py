#!/usr/bin/env python

import os
import subprocess
from distutils.core import setup

# Get the git version
try:
    git_version = subprocess.check_output('git describe --abbrev=6 --dirty --always', shell=True).strip()
except:
    git_version = '?'

# Store it in the unversioned file _version.
# This could be handled using a Makefile target
fh = open('atooms/_version.py', 'w')
fh.write('__version__ = "%s"' % git_version)
fh.close()

setup(name='atooms',
      version=git_version,
      description='Atomistic object-oriented modeling and simulations in Python',
      author='Daniele Coslovich',
      author_email='daniele.coslovich@univ-montp2.fr',
      url='some.site.fr',
      packages=[],
      # 'atooms', 'atooms.adapters', 'atooms.system',
      # 'atooms.interaction', 'atooms.potential', 
      # 'atooms.postprocessing'],
      scripts=['bin/trj.py']
      #scripts=glob.glob(os.path.join('bin', '*.py'))
     )
