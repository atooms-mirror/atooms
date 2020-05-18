#!/usr/bin/env python

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

# Get the long description from README.md and try to convert it to
# reST. Adapted from https://bons.ai/blog/markdown-for-pypi
try:
    from pypandoc import convert
    readme = convert('README.md', 'rst')
except (ImportError, OSError):
    try:
        readme = open('README.md', 'r').read()
    except:
        readme = ''

with open('atooms/core/_version.py') as f:
    exec(f.read())

setup(name='atooms',
      version=__version__,
      description='A framework for simulations of interacting particles',
      long_description=readme,
      author='Daniele Coslovich',
      author_email='daniele.coslovich@umontpellier.fr',
      url='https://gitlab.info-ufr.univ-montp2.fr/atooms/atooms',
      packages=['atooms', 'atooms/backends', 'atooms/backends/f90', 'atooms/core', 'atooms/interaction', 
                'atooms/plugins', 'atooms/simulation', 'atooms/system', 
                'atooms/trajectory'],
      scripts=['bin/trj.py'],
      long_description_content_type="text/markdown",
      install_requires=['numpy'],
      package_data = {'atooms/backends/f90': ['*.f90']},
      license='GPLv3',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Physics',
      ]
)
