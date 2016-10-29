#!/usr/bin/env python

from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

setup(name='atooms',
      version='0.1',
      description='A framework for classical particle-based simulations',
      long_description=readme,
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms',
      packages=['atooms',
                'atooms.backends', 
                'atooms.interaction',
                'atooms.potential',
                'atooms.simulation',
                'atooms.system', 
                'atooms.trajectory'],
      license='GPLv3',
      classifiers={
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
      }
)
