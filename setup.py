#!/usr/bin/env python

import os
import subprocess
from distutils.core import setup

# Get the git version
try:
    git_version = subprocess.check_output('git describe --abbrev=6 --dirty --always', shell=True).strip()
    git_date = subprocess.check_output('git show -s --format=%ci %s' % git_version, shell=True).strip()
except:
    git_version = '?'
    git_date = '?'

# Store it in the unversioned file _version.
# This could be handled using a Makefile target
fh = open('atooms/_version.py', 'w')
fh.write('__version__ = "%s"' % git_version)
fh.write('__date__ = "%s"' % git_date)
fh.close()

setup(name='atooms',
      version=git_version,
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
