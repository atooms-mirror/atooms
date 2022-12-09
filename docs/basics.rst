


Basics
------

Atooms provides a high-level interface to the main objects of particle-based simulations. It mostly focuses on classical molecular dynamics and Monte Carlo simulations, but it is not limited to that. It can be used to simulate and analyze lattice models such as TASEP or kinetically constrained models.

We will start by having a look at the basic objects of particle-based simulations and how to store them on a file.

Particles' properties
~~~~~~~~~~~~~~~~~~~~~

Particles' positions are stored as numpy arrays, but we can pass a simple list with x, y, z coordinates when we create them

.. code:: python

    from atooms.system.particle import Particle
    particle = Particle(position=[1.0, 0.0, 0.0])
    print(particle.position, type(particle.position))

::

    Python 3.8.10 (default, Jun 22 2022, 20:18:18) 
    [GCC 9.4.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    [1. 0. 0.] <class 'numpy.ndarray'>

Particles can live in an arbitrary number of spatial dimensions

.. code:: python

    particle = Particle(position=[1.0, 0.0, 0.0, 0.0, 0.0])
    print(len(particle.position))

::

    5

By default, particles have a few more properties such as velocity, chemical species, mass and radius. They can all be altered at will or even set to None.

.. code:: python

    import numpy
    particle = Particle(position=[1.0, 0.0, 0.0], velocity=[1.0, 0.0, 0.0])
    particle.species = 'Na'
    particle.position += numpy.array([0.0, 1.0, 1.0])
    particle.velocity *= 2
    particle.radius = None  # point particles have no radius
    print(particle)

::

    Particle(species=Na, mass=1.0, position=[1. 1. 1.], velocity=[2. 0. 0.], radius=None)

You may want to add physical properties to particles, like charge or whatever. Of course, in python you can do it very easily

.. code:: python

    particle.charge = -1.0

This won't break anything!

Dealing with velocities
~~~~~~~~~~~~~~~~~~~~~~~

You may not need velocities at all (for instance because you are working with Monte Carlo simulations) but if you do, atooms provides a few useful methods and functions. For instance, you can assign velocity from a Maxwell-Boltzmann distribution at a temperature T.

.. code:: python

    particle = [Particle() for i in range(1000)]
    for p in particle:
        p.maxwellian(T=1.0)
    ekin = sum([p.kinetic_energy for p in particle])
    ndim = 3
    ndof = len(particle) * ndim
    T = 2.0 / ndof * ekin
    print(T)

::

    1.0109937775682067

Doing so will leave a non-zero total momentum, but we can fix it (note that all masses are equal)

.. code:: python

    from atooms.system.particle import fix_total_momentum, cm_velocity
    print(cm_velocity(particle))
    fix_total_momentum(particle)
    print(cm_velocity(particle))

::

    [-0.02489158 -0.05049646  0.00014979]
    [-1.15463195e-17  4.34097203e-17  6.99440506e-18]

Boundary conditions
~~~~~~~~~~~~~~~~~~~

To avoid major finite size effects, we enclose particles in a cell with periodic boundary conditions. By convention, the cell origin is at the origin of the reference frame.

.. code:: python

    from atooms.system.cell import Cell
    L = 2.0
    cell = Cell(side=[L, L, L])
    print(cell.side, cell.volume)

::

    [2. 2. 2.] 8.0

Atooms provides means to fold particles back in the "central" simulation cell, i.e. the one centered at the origin at the reference frame. For simplicity, let us work with particles in 1d.

.. code:: python

    cell = Cell(side=1.0)
    particle = Particle(position=2.0)  # particle outside the central cell
    particle.fold(cell)
    print(particle.position)

::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/tmp/python-FjISQM", line 3, in <module>
        particle.fold(cell)
      File "/home/coslo/envs/dev/lib/python3.8/site-packages/atooms/system/particle.py", line 95, in fold
        self.position[:] = _periodic_vector_unfolded(self.position, cell.side)
    IndexError: too many indices for array: array is 0-dimensional, but 1 were indexed

The particle is now folded back at the origin.

A related method returns the nearest periodic image of a given particle with respect to another particle

.. code:: python

    particle_1 = Particle(position=-0.45)
    particle_2 = Particle(position=+0.45)
    image = particle_1.nearest_image(particle_2, cell, copy=True)
    print(image)

::

    Particle(species=A, mass=1.0, position=0.55, velocity=[0. 0. 0.], radius=0.5)

The System object
~~~~~~~~~~~~~~~~~

Objects like particles and the simulation cell can be gathered in an instance of a god-like class called System. The system contains all the relevant physical objects of your simulation. Reservoirs like thermostats, barostats and particle reservoirs can be added as well. These objects are placeholders for thermodynamic state variables like temperature, pressure or chemical potential. Any class meant to describe the interaction between particles also belongs to the system.

Let us build a system with a few particles in a cell and use the system methods to modify the system density and temperature. Note that density and temperature are python properties and thus modify the attributes of particles and cell under the hoods using the ``set_density`` and ``set_temperature`` methods respectively

.. code:: python

    from atooms.system import System
    system = System(particle=[Particle() for i in range(100)],
    		cell=Cell([10.0, 10.0, 10.0]))
    system.density = 1.2  # equivalent to system.set_density(1.2)
    system.temperature = 1.5  # equivalent to system.set_temperature(1.2)
    print(system.density, system.temperature)

::

    1.1999999999999997 1.5000000000000004

Note that the system temperature is the kinetic one and need not coincide with the one of the thermostat.

.. code:: python

    from atooms.system import Thermostat
    system.thermostat = Thermostat(temperature=1.0)
    system.temperature = 1.5  # equivalent to system.set_temperature(1.2)
    print(system.temperature, system.thermostat.temperature)

::

    1.5000000000000002 1.0

Interaction and backends
~~~~~~~~~~~~~~~~~~~~~~~~

Classical particles interact with each other via a potential :math:`u(\{r_i\})`, where :math:`\{r_i\}` is the set of particles' coordinates. Atooms relies on third-party efficient **backends** written in C, Fortran or CUDA to actually compute the interaction between the particles. Here we will use the LAMMPS backend, see Molecular dynamics ith LAMMPS for further details. It accepts a string variable that defines the interaction potential using the LAMMPS syntax, see `https://lammps.sandia.gov/doc/pair_style.html <https://lammps.sandia.gov/doc/pair_style.html>`_, and stores a reference to the system object of which we want to compute the energy.

As proof of principle, we compute the interaction energy between two Lennard-Jones particles

.. code:: python

    from atooms.system import System, Particle, Cell
    from atooms.backends.lammps import LAMMPS

    x = 1.122  # Minimum of the potential
    system = System(particle=[Particle(position=[0.0, 0.0, 0.0]),
    			  Particle(position=[x, 0.0, 0.0])],
    		cell=Cell([10.0, 10.0, 10.0]))
    cmd = """
    pair_style      lj/cut 2.5
    pair_coeff      1 1 1.0 1.0  2.5
    """
    # The backend will add an interaction to the system
    backend = LAMMPS(system, cmd)

    # Compute and get the potential energy
    # The cache option allows to get the potential energy without recalculating it
    print(system.potential_energy(), system.potential_energy(cache=True))

::

    Traceback (most recent call last):
      File "/home/coslo/envs/dev/lib/python3.8/site-packages/atooms/backends/lammps.py", line 48, in _get_lammps_version
        _ = subprocess.check_output(cmd, shell=True,
      File "/usr/lib/python3.8/subprocess.py", line 415, in check_output
        return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
      File "/usr/lib/python3.8/subprocess.py", line 516, in run
        raise CalledProcessError(retcode, process.args,
    subprocess.CalledProcessError: Command 'echo | mpirun lammps' returned non-zero exit status 134.

    During handling of the above exception, another exception occurred:

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/tmp/python-R17fkG", line 13, in <module>
        backend = LAMMPS(system, cmd)
      File "/home/coslo/envs/dev/lib/python3.8/site-packages/atooms/backends/lammps.py", line 158, in __init__
        self.version = _get_lammps_version()
      File "/home/coslo/envs/dev/lib/python3.8/site-packages/atooms/backends/lammps.py", line 52, in _get_lammps_version
        raise ImportError('lammps not installed (command is {})'.format(lammps_command))
    ImportError: lammps not installed (command is lammps)

The energy and forces are stored in ``system.interaction.energy`` and ``system.interaction.forces``.

Trajectory files
~~~~~~~~~~~~~~~~

To write the state of the system to a file, we use a ``Trajectory`` class. Trajectories are composed of multiple frames, each one holding the state of the system at a given step during the simulation. We use a basic xyz format to write the state of the system and then parse the trajectory file we produced to see how it looks like.

.. code:: python

    from atooms.trajectory import TrajectoryXYZ

    system = System(particle=[Particle() for i in range(4)],
                    cell=Cell([10.0, 10.0, 10.0]))

    # Open the trajectory in write mode and write the state of the system
    # at step 0
    with TrajectoryXYZ('test.xyz', 'w') as th:
        th.write(system, step=0)

    # Read the xyz file back as plain text
    with open('test.xyz') as fh:
        print(fh.read())

::

    4
    step:0 columns:species,position dt:1 cell:10.0,10.0,10.0 
    A 0.000000 0.000000 0.000000
    A 0.000000 0.000000 0.000000
    A 0.000000 0.000000 0.000000
    A 0.000000 0.000000 0.000000

Note that trajectories are file-like objects: they must be opened and closed, preferably using the ``with`` syntax.

Of course, we can write multiple frames by calling ``write()`` repeatedly.

.. code:: python

    with TrajectoryXYZ('test.xyz', 'w') as th:
        for i in range(3):
            th.write(system, step=i*10)

To get the system back we read the trajectory. Trajectories support iteration and indexing, just like lists.

.. code:: python

    with TrajectoryXYZ('test.xyz') as th:
        # First frame
        system = th[0]
        print(system.particle[0].position, system.cell.side)

        # Last frame
        system = th[-1]
        print(system.particle[0].position, system.cell.side)

        # Iterate over all frames
        for i, system in enumerate(th):
            print(th.steps[i], system.particle[0].position)

::

    [0. 0. 0.] [10. 10. 10.]
    [0. 0. 0.] [10. 10. 10.]
    0 [0. 0. 0.]
    10 [0. 0. 0.]
    20 [0. 0. 0.]

Particles on a lattice
~~~~~~~~~~~~~~~~~~~~~~

Suppose we want to simulate a system where particles can only be located at discrete sites, say a one-dimensional lattice or perhaps a network with a complex topology. Particle positions can then be described as plain integers, holding the index of the site on which a particle is located. We create such a system and then write it to a file in xyz format

.. code:: python

    import numpy
    from atooms.system import System, Particle

    # Build model system with integer coordinates
    particle = [Particle() for i in range(3)]
    particle[0].position = 0
    particle[1].position = 1
    particle[2].position = 2
    system = System(particle=particle)

    # Write xyz trajectory
    from atooms.trajectory import TrajectoryXYZ
    with TrajectoryXYZ('test.xyz', 'w') as th:
        th.write(system, 0)

    # Read the xyz file back as plain text
    with open('test.xyz') as fh:
        print(fh.read())

::

    3
    step:0 columns:species,position dt:1 
    A 0
    A 1
    A 2

Everything went fine. However, we have to tweak things a bit when reading the particles back, to avoid positions being transformed to arrays of floats instead of integers. This can be done with the help of a callback that transforms the system accordingly as we read the trajectory.

.. code:: python

    # Read file as an xyz trajectory 
    with TrajectoryXYZ('test.xyz') as th:
        # We add a callback to read positions as simple integers
        # Otherwise they are read as numpy arrays of floats.
        def modify(system):
            for p in system.particle:
                p.position = int(p.position[0])
                p.velocity = None
                p.radius = None
            return system
        th.add_callback(modify)

        for p in th[0].particle:
            print(p)

::

    Particle(species=A, mass=1.0, position=0, velocity=None, radius=None)
    Particle(species=A, mass=1.0, position=1, velocity=None, radius=None)
    Particle(species=A, mass=1.0, position=2, velocity=None, radius=None)

Our particles have now integer coordinates. Note that, on passing, we have set to None velocities and radii as they are not relevant in this case.
