#!/usr/bin/env python

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

with open('atooms/core/_version.py') as f:
    exec(f.read())

setup(name='atooms',
      version=__version__,
      description='A framework for classical simulations of interacting particles',
      long_description=readme,
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms',
      packages=['atooms', 'atooms/core', 'atooms/interaction', 'atooms/plugins', 
                'atooms/simulation', 'atooms/system', 'atooms/trajectory', 
                'atooms/utils'],
      scripts=['bin/trj.py', 'bin/md.py'],
      install_requires=['numpy'],
      license='GPLv3',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
      ]
)
