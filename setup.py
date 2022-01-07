#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md') as f:
    readme = f.read()

with open('atooms/core/_version.py') as f:
    exec(f.read())

setup(name='atooms',
      version=__version__,
      description='A framework for simulations of interacting particles',
      long_description=readme,
      long_description_content_type="text/markdown",
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://framagit.org/atooms/atooms',
      packages=['atooms', 'atooms/backends', 'atooms/core',
                'atooms/plugins', 'atooms/simulation', 'atooms/system', 
                'atooms/trajectory'],
      scripts=['bin/trj.py'],
      install_requires=['numpy'],
      license='GPLv3',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.9',
          'Topic :: Scientific/Engineering :: Physics',
      ]
)
