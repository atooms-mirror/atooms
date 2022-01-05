# Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). This file only reports changes that increase major and minor versions, as well as deprecations.

## 3.5.0 - 2022/01/05

[Full diff](https://framagit.org/atooms/atooms/-/compare/3.4.3...3.5.0)

## Added
- Add system.cm() to compute general attribute of CM

## Changed
- Remove `particle` from attributes in comment header of xyz files. This is backward compatible.

## 3.4.0 - 2021/11/18

[Full diff](https://framagit.org/atooms/atooms/-/compare/3.3.4...3.4.0)

## Added
- Refactor and improve simulation observers
  - Refactor `write_thermo()` and `write_config()` as stateless functions
  - Add optional `trajectory_class` parameter to `write_trajectory()` to choose the trajectory format
  - Add optional `trajectory` parameter to `write_trajectory()` to write to an existing `Trajectory`
  - Add `store()` simulation observer to store properties in a dict during the simulation
  - Add option `variables` to `write_trajectory()` simulation observe
  - Support more general attributes writing in `write()` via the `what` variable
  - Expect deprecations for `write_config()` (use `write_trajectory` instead)
  - Expect deprecations for `write_to_ram()` (use `write_trajectory` passing a `TrajectoryRAM` instance instead)
- Add default value for `Simulation.trajectory_class` as `TrajectoryXYZ`
- Add `coordinate` and `momentum` arrays to reservoirs to enable in-place modification (f2py) and chains
- Reservoir masses are arrays too
	 
## 3.3.0 - 2021/11/02

[Full diff](https://framagit.org/atooms/atooms/-/compare/3.2.0...3.3.0)

## Added

- Add `Wall` class and optional instances as `System.wall`
- Add `System.species_layout` property to show and change the chemical species layout (A, C, F)
- Add `InteractionBase` as the base for actual Interaction subclasses
- Refactor `Interaction` to handle multiple interaction terms via the `term` list variable

## 3.2.0 - 2021/10/31

[Full diff](https://framagit.org/atooms/atooms/-/compare/3.1.0...3.2.0)

### Added

- Add syntax to optionally specify the data type of Interaction.variables as <property[:dtype]>. This allows dumping arrays with the correct data type.

## 3.1.0 - 2021/10/30

[Full diff](https://framagit.org/atooms/atooms/-/compare/3.0.0...3.1.0)

### Added

- Add `TrajectoryFactory.register_callback()`
- Argument `what` in `System.compute_interaction()` is now optional

## 3.0.0 - 2021/10/03

[Full diff](https://framagit.org/atooms/atooms/-/compare/2.8.1...3.0.0)
	
### Changed
- Interaction class now belongs to system.interaction module
- Change the interface of Interaction.compute() method so as to accept System arbitrary attributes. Thanks to this change, Interaction now behaves as a strategy pattern for System
- Use new Trajectory.copy() method to convert trajectory, instead of trajectory.utils.convert()
- Rename inappropriate Trajectory.trajectory attribute in Trajectory subclasses. The new convention is to use _file as private attribute for the file object
- Rename index attribute to frame in TrajectoryBase.read_system()
- Change options order in Particle constructor
- Change system.distinct_species() method as a property
- RUMD backend requires version 3
- Rename forcefield option as potentials in RUMD backend
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
- Remove forcefield_file option from RUMD backend, pass RUMD Potentials instances instead
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


