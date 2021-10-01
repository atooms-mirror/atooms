# Changelog

The format of this file is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

[Full changelog](https://framagit.org/atooms/atooms/-/compare/2.8.1...master)
	
### Changed
- Interaction class now belongs to system.interaction module
- Change the interface of Interaction.compute() method so as to accept System arbitrary attributes. Thanks to this change, Interaction now behaves as a strategy pattern for System
- Use new Trajectory.copy() method to convert trajectory, instead of trajectory.utils.convert()
- Rename inappropriate Trajectory.trajectory attribute in Trajectory subclasses. The new convention is to use _file as private attribute for the file object
- Rename index attribute to frame in TrajectoryBase.read_system()
- Change options order in Particle constructor
- Change system.distinct_species() method as a property
- Merge autopep8 and flake8 Makefile targets into pep8

### Removed
- Remove Interaction subpackage, including potential tabulation stuff, which will be moved to a specific backend.
- Remove atooms.backends.f90, it will be added to the atooms-models package
- Remove Trajectory.self_callbacks and Trajectory.class_callbacks
- Remove atooms.simulation.umbrella
- Remove system.particle.show() function, because it belongs to System
- Remove system.particle.show_*() functions, because they belong to visualize
- Remove trajectory.utils.convert() and trajectory.utils.modify_fields()
- Remove unfolded parameter from TrajectoryHDF5.read_sample()
- Remove support for cell values as last line in xyz trajectory format
- Remove support for matrix attribute in System
- Remove support for list of attributes passed to System.dump(), use dict compreheshion instead to collect dumped attributes
- Remove report() methods throughout, rely on __str__() instead
- Remove user and develop targets from Makefile
- Remove tox.ini

### Deprecated
- Trajectory.write_sample() in favor of write_system()
- Trajectory.read_sample() in favor of read_system()
- Trajectory.fields in favor of Trajectory.variables

### Added
- Add Trajectory.variables to store System attributes that change along the trajectory
- Add core.utils.canonicalize() as standalone function
- Add contribution guidelines
- Add setup.cfg configuration file


