# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Interaction class is now in system.interaction module
- Change the interface of Interaction compute() method to accept System arbitrary attributes
- Refactor Interaction as a strategy in System
- Use Trajectory.copy() to convert trajectories instead of trajectory.utils.convert()
- Privatize and rename inappropriate Trajectory.trajectory attribute in Trajectory subclasses: the new convention is to use _file as private attribute for the file object
- Rename index attribute to frame in TrajectoryBase.read_system()
- Change options order in Particle.__init__()
- Change system.distinct_species() method into a property
- Merge autopep8 and flake8 Makefile targets into pep8

### Removed
- Remove Interaction subpackage
- Remove atooms.backends.f90, it will be added to the atooms-models package
- Remove self_callbacks and class_callbacks from Trajectory
- Remove atooms.simulation.umbrella`
- Remove show() and show_*() from system.particle, they belong to System
- Remove convert() and modify_fields() from trajectory.utils
- Remove unfolded parameter from TrajectoryHDF5.read_sample()
- Remove support for cell values as last line in xyz trajectory format
- Remove support for matrix attribute in System
- Remove support for list of attributes passed to System.dump(): use dict compreheshion instead to collect dumped attributes
- Remove report() methods throughout: rely on __str__() instead
- Remove user and develop targets from Makefile
- Remove tox.ini

### Deprecated
- Deprecate Trajectory.write_sample() in favor of write_system()
- Deprecate Trajectory.read_sample() in favor of read_system()
- Deprecate Trajectory.fields in favor of Trajectory.variables

### Added
- Add Trajectory.variables to store System attributes that change along the trajectory
- Add core.utils.canonicalize() as standalone function
- Add contribution guidelines
- Add setup.cfg configuration file

### Fixed
- Small fixes
- Code refactoring


