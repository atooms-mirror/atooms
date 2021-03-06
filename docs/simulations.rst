


Simulations
-----------

atooms provides a generic interface that abstracts out most of the common tasks of particle-based simulations. The actual simulation is performed by a simulation backend, which exposes a minimal but consistent interface. This enables one to develop complex simulation frameworks (e.g., [parallel tempering](`https://framagit.org/atooms/parallel_tempering <https://framagit.org/atooms/parallel_tempering>`_)) that are essentially decoupled from the underlying simulation code.

A **Simulation** is a high-level class that encapsulates some common tasks and provides a consistent interface to the user, while **backend** classes actually make the system evolve. Here, we implement a minimal backend to run a simulation.

At a very minimum, a backend is a class that provides 

- a **system** instance variable, which should (mostly) behave like ``atooms.system.System``.

- a **run()** method, which evolves the system for a prescribed number of steps (passed as argument)

Optionally, the backend may hold a reference to a trajectory class, which can be used to checkpoint the simulation or to write configurations to a file. This is however not required in a first stage.

A minimal simulation backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We set up a bare-bones simulation backend building on the native System class

.. code:: python

    from atooms.system import System

    class BareBonesBackend(object):
    
        def __init__(self):
            self.system = System()

        def run(self, steps):
            for i in range(steps):
                pass

    # The backend is created and wrapped by a simulation object.
    # Here we first call the run() method then run_until()
    from atooms.simulation import Simulation
    backend = BareBonesBackend()
    simulation = Simulation(backend)
    simulation.run(10)
    simulation.run_until(30)
    assert simulation.current_step == 30

    # This time we call run() multiple times 
    simulation = Simulation(backend)
    simulation.run(10)
    simulation.run(20)
    assert simulation.current_step == 30  

    # Increase verbosity to see a meaningful log
    from atooms.core.utils import setup_logging
    setup_logging(level=20)
    simulation = Simulation(backend)
    simulation.run(10)  

::

    # 
    # atooms simulation via <__main__.BareBonesBackend object at 0x7ff54d0527f0>
    # 
    # version: 1.9.1+1.5.0-132-gfe9bc7-dirty (2019-04-12)
    # atooms version: 1.9.1+1.5.0-132-gfe9bc7-dirty (2019-04-12)
    # simulation started on: 2019-05-17 at 17:36
    # output path: None
    # backend: <__main__.BareBonesBackend object at 0x7ff54d0527f0>
    # 
    # target target_steps: 10
    # 
    # 
    # starting at step: 0
    # 
    # simulation ended successfully: reached target steps 10
    # 
    # final steps: 10
    # final rmsd: 0.00
    # wall time [s]: 0.00
    # average TSP [s/step/particle]: nan
    # simulation ended on: 2019-05-17 at 17:36

Simple random walk
~~~~~~~~~~~~~~~~~~

We implement a simple random walk in 3d. This requires adding code to the backend ``run()`` method to actually move the particles around.

We start by building an empty system. Then we add a few particles and place them at random in a cube. Finally, we write a backend that displaces each particle randomly over a cube of prescribed side.

.. code:: python

    import numpy
    from atooms.system import System

    # There are no particles at the beginning
    system = System()
    assert len(system.particle) == 0

    # Add particles
    from atooms.system.particle import Particle
    from random import random
    L = 10
    for i in range(1000):
        p = Particle(position=[L * random(), L * random(), L * random()])
        system.particle.append(p)

    class RandomWalk(object):

        def __init__(self, system, delta=1.0):
            self.system = system
            self.delta = delta

        def run(self, steps):
            for i in range(steps):
                for p in self.system.particle:
                    dr = numpy.array([random()-0.5, random()-0.5, random()-0.5])
                    dr *= self.delta
                    p.position += dr

The Simulation class provides a callback mechanism to allow execution of arbitrary code during the simulation. This can be used to write logs or particle configurations to file, or to perform on-the-fly calculations of the system properties. Callbacks are plain function that accept the simulation object as first argument. They are called at prescribed intervals during the simulation.

Here we measure the mean square displacement (MSD) of the particles to make sure that the system displays a regular diffusive behavior :math:`MSD \sim t`

.. code:: python

    from atooms.simulation import Simulation
    simulation = Simulation(RandomWalk(system))

    # We add a callback that computes the MSD every 10 steps
    # We store the result in a dictionary passed to the callback
    msd_db = {}
    def cbk(sim, initial_position, db):
        msd = 0.0
        for i, p in enumerate(sim.system.particle):
            dr = p.position - initial_position[i]
            msd += numpy.sum(dr**2)
        msd /= len(sim.system.particle)
        db[sim.current_step] = msd

    # We will execute the callback every 10 steps
    simulation.add(cbk, 10, initial_position=[p.position.copy() for p in
                                              system.particle], db=msd_db)
    simulation.run(50)

    # The MSD should increase linearly with time
    time = sorted(msd_db.keys())
    msd = [msd_db[t] for t in time]

    print(time, msd)
    import matplotlib.pyplot as plt
    plt.cla()
    plt.plot(time, msd, '-o')
    plt.xlabel("t")
    plt.ylabel("MSD")
    plt.savefig('msd.png')

The MSD as a function of time should look linear.
.. image:: msd.png

Molecular dynamics with LAMMPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Atooms provides a simulation backend for ``LAMMPS``, an efficient and feature-rich molecular dynamics simulation package.
The backend accepts a string variable containing regular LAMMPS commands and initial configuration to start the simulation. The latter can be provided in any of the following forms:

- a System object

- a Trajectory object

- the path to an xyz trajectory

In the last two cases, the last configuration will be used to start the simulation. 

Here we we use the first configuration of an existing trajectory for a Lennard-Jones fluid

.. code:: python

    import atooms.trajectory as trj
    from atooms.backends.lammps import LAMMPS

    import os
    system = trj.TrajectoryXYZ('../../data/lj_N1000_rho1.0.xyz')[0]
    cmd = """
    pair_style      lj/cut 2.5
    pair_coeff      1 1 1.0 1.0  2.5
    neighbor        0.3 bin
    neigh_modify    check yes
    timestep        0.002
    """
    backend = LAMMPS(system, cmd)

We now wrap the backend in a simulation instance. This way we can rely on atooms to write thermodynamic data and configurations to disk during the simulation: we just add the ``write_config()`` and ``write_thermo()`` functions as observers to the simulations.
You can add your own functions as observers to perform arbitrary manipulations on the system during the simulation. Keep in mind that calling these functions causes some overhead, so avoid calling them at too short intervals.

.. code:: python

    from atooms.simulation import Simulation
    from atooms.system import Thermostat
    from atooms.simulation.observers import write_thermo, write_config

    # We create the simulation instance and set the output path
    sim = Simulation(backend, output_path='lammps.xyz')
    # Just store a reference to the trajectory class you want to use
    sim.trajectory_class = trj.TrajectoryXYZ
    # Write configurations every 500 steps in xyz format
    sim.add(write_config, 500)
    # Write thermodynamic properties every 500 steps
    sim.add(write_thermo, 500)

We add a thermostat to keep the system temperature at T=2.0 and run the simulations for 10000 steps.

.. code:: python

    backend.system.thermostat = Thermostat(temperature=2.0, relaxation_time=0.1)
    sim.run(10000)

Note that we use atooms ``Thermostat`` object here: the backend will take care of adding appropriate commands to the LAMMPS script.

We have a quick look at the kinetic temperature as function of time to make sure the thermostat is working

.. code:: gnuplot

    set xl 'Steps'
    set yl 'Temperature'
    set border 3
    set xtics nomirror
    set ytics nomirror
    plot 'lammps.xyz.thermo' u 1:2 noti w lp lc rgb 'red' pt 7, 2 noti lc rgb 'black'

.. image:: lammps.png

We can use the `postprocessing <https://gitlab.info-ufr.univ-montp2.fr/atooms/postprocessing/>`_ package to compute the radial distribution function

.. code:: python

    from atooms.postprocessing import api
    api.gr('lammps.xyz')

.. code:: gnuplot

    set xl 'r'
    set yl 'g(r)'
    set border 3
    set xtics nomirror
    set ytics nomirror
    plot 'lammps.xyz.pp.gr' u 1:2 noti w lp lc rgb 'red' pt 7

.. image:: lammps_gr.png

Molecular dynamics simulation with RUMD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we pick the last frame of the trajectory, change the density of the system to unity and write this new configuration to a trajectory format suitable for the `RUMD <http://rumd.org>`_ simulation package

.. code:: python

    with Trajectory('input.xyz') as trajectory:
        system = trajectory[-1]
        system.density = 1.0
        print('New density:', len(system.particle) / system.cell.volume)

    from atooms.trajectory import TrajectoryRUMD
    with TrajectoryRUMD('rescaled.xyz.gz', 'w') as trajectory:
        trajectory.write(system)

Now we run a short molecular dynamics simulation with the ``RUMD`` backend, using a Lennard-Jones potential:

.. code:: python

    import rumd
    from atooms.backends.rumd import RUMD
    from atooms.simulation import Simulation

    potential = rumd.Pot_LJ_12_6(cutoff_method=rumd.ShiftedPotential)
    potential.SetParams(i=0, j=0, Epsilon=1.0, Sigma=1.0, Rcut=2.5)
    backend = RUMD('rescaled.xyz.gz', [potential], integrator='nve'
    sim = Simulation(backend)
    sim.run(1000)
    print('Final temperature and density:', sim.system.temperature, sim.system.density)

Energy minimization with LAMMPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to minimize the energy of a system to determine its so-called inherent structure using LAMMPS as a backend. To achieve this, atooms defines an ``Optimization`` class, which behaves mostly as ``Simulation`` except that it stops when the mean square total force


.. math::

    W=\frac{1}{N}\sum_i |f_i|^2


is lower than a given ``tolerance``.

.. code:: python

    from atooms.trajectory import TrajectoryXYZ
    from atooms.optimization import Optimization
    from atooms.backends.lammps import EnergyMinimization
    cmd = """
    pair_style      lj/cut 2.5
    pair_modify     shift yes
    pair_coeff      1 1 1.0 1.0 2.5
    """
    system = TrajectoryXYZ('../../data/lj_N256_rho1.0.xyz')[0]
    bck = EnergyMinimization(system, cmd)
    opt = Optimization(bck, tolerance=1e-10)
    opt.run()

We check that :math:`W` is lower than the requested tolerance

.. code:: python

    e_final = system.potential_energy(per_particle=True)
    w_final = system.force_norm_square(per_particle=True)
    print('Energy={}, mean square force={:.2g}'.format(e_final, w_final))

::

    Energy=-6.8030584, mean square force=3.6e-11
