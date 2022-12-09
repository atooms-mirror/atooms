


Trajectories
------------

Trajectory formats
~~~~~~~~~~~~~~~~~~

``atooms`` supports several trajectory formats, most of them in read and write mode

.. code:: python

    from atooms.trajectory import formats

    print(formats())

::

    Available trajectory formats:
    - dynamo       : [R ] DynamO trajectory format (https://www.dynamomd.com/index.php/tutorial3)
    - exyz         : [RW] Extended XYZ layout (https://github.com/libAtoms/extxyz)
    - folderlammps : [R ] Multi-file layout LAMMPS format.
    - gsd          : [RW] Glotzer group's binary GSD format for HOOMD (https://glotzerlab.engin.umich.edu/hoomd-blue/)
    - hdf5         : [RW] In-house trajectory layout in HDF5 format.
    - hoomd        : [RW] HOOMD format
    - lammps       : [RW] LAMMPS format (https://docs.lammps.org/dump.html)
    - neighbors    : [RW] Neighbors trajectory format
    - pdb          : [RW] PDB format (https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format))
    - ram          : [RW] Store trajectory in RAM
    - rumd         : [RW] RUMD trajectory format (https://rumd.org)
    - simplexyz    : [RW] Simple implementation of the xyz layout (https://en.wikipedia.org/wiki/XYZ_file_format)
    - superrumd    : [R ] SuperTrajectory for RUMD format
    - xyz          : [RW] XYZ format with metadata support (https://en.wikipedia.org/wiki/XYZ_file_format)

Custom trajectory formats
~~~~~~~~~~~~~~~~~~~~~~~~~

It is easy to add new trajectory formats by subclassing existing trajectory classes. To make these new classes accessible also to ``trj.py``, create a package called ``atooms_plugins`` and add your trajectory modules there. Suppose you wrote a custom trajectory class ``TrajectoryABC`` in ``atooms_plugins/test.py`` (the last path is relative to the current directory). You can now convert an existing xyz trajectory to your custom format:

.. code:: sh

    trj.py convert output.xyz output.abc

Remember to add an empty ``__init__.py`` file at the root of ``atooms_plugins``. Actually, the ``atooms_plugins`` package can be put anywhere in your ``PYTHONPATH``.

Custom trajectory output
~~~~~~~~~~~~~~~~~~~~~~~~

We can customize the format of trajectory files using the ``fields`` variable. It contains a list of the particle properties to be written to the trajectory. For this simple example we use again the xyz trajectory format.

We add a ``charge`` property to each particle and then instruct the trajectory to write it along with the position

.. code:: python

    from atooms.system import System, Cell, Particle
    system = System(particle=[Particle() for i in range(3)],
    		cell=Cell([10.0, 10.0, 10.0]))

    for p in system.particle:
        p.charge = -1.0

    with TrajectoryXYZ('test.xyz', 'w', fields=['position', 'charge']) as th:
        th.write(system, step=0)

    with open('test.xyz') as fh:
        print(fh.read())

::

    3
    step:0 columns:position,charge dt:1 cell:10.0,10.0,10.0 
    0.000000 0.000000 0.000000 -1.0
    0.000000 0.000000 0.000000 -1.0
    0.000000 0.000000 0.000000 -1.0

The ``fields`` list can contain any particle property, even those defined dynamically at run time, such as the ``charge`` variable above which is not a predefined particle property!. When reading back the trajectory, the ``charge`` property is automatically recognized and added to the particle. 

.. code:: python

    with TrajectoryXYZ('test.xyz') as th:
      system = th[0]
      print(system.particle[0].charge)

::

    -1.0

Conversion between trajectory formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Atooms provides means to convert between trajectory various formats. At a very basic level, this requires opening the original trajectory for reading and the new one for writing using the desired trajectory class. Here we convert an xyz trajectory in a format suitable for the LAMMPS package

.. code:: python

    from atooms.trajectory import TrajectoryLAMMPS
    with TrajectoryXYZ('test.xyz') as th_inp,\
         TrajectoryLAMMPS('test.lammps', 'w') as th_out:
        for i, system in enumerate(th_inp):
            th_out.write(system, th_inp.steps[i])

The ``convert()`` function wraps the conversion in a more convenient interface

.. code:: python

    from atooms.trajectory import convert
    convert(TrajectoryXYZ('test.xyz'), TrajectoryLAMMPS, 'test.lammps')

There are several optional parameters that allows to customize the trajectory output, see the function signature for more details.

Finally, the ``trj.py`` script installed by atooms allows to quickly convert trajectories on the command-line, which is actually the most frequent use case

.. code:: sh

    trj.py convert -i xyz -o lammps test.xyz test.lammps

Although the script will do its best to guess the appropriate trajectory formats, it is best to provide the input and output trajectory formats via the ``-i`` and ``-o`` flags explicitly.

Add and modify trajectory properties on the fly with callbacks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"Callbacks" are functions used to modify the properties of a trajectory on the fly. They accept a ``System`` instance as first positional argument, along with optional extra positional and keyword arguments, and return a modified ``System``.

As an example, suppose your trajectory did not provide any information about the cell side. You can add the information dynamically to all ``System`` objects read from the trajectory using the following callback

.. code:: python

    from atooms.system import Cell
    def fix_missing_cell(system, side):
        system.cell = Cell(side)
        return system

Then we add the callback to the trajectory and provide the cell side (here L=10 along each dimensions) as argument. Reading the trajectory is then done as usual.

.. code:: python

    from atooms.trajectory import TrajectoryXYZ
    with TrajectoryXYZ('test.xyz') as th:
        th.add_callback(fix_missing_cell, [10., 10., 10.])
        for system in th:
            print(system.cell.side)

::

    [10. 10. 10.]
    [10. 10. 10.]

Extend trajectory classes
~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you have a trajectory that looks almost like xyz, but differs in some way. You may want to customize the xyz trajectory format, so that your code can process the trajectory without manual intervention.

For instance, your xyz file is ``test.xyz`` but the cell side information is stored in a separate file ``test.xyz.cell``. We can proceed as before

.. code:: python

    from atooms.system import Cell

    file_inp = 'test.xyz'
    with open(file_inp + '.cell') as fh:
        # Assume the cell file contains a string Lx Ly Lz
        # where Lx, Ly, Lz are the sides of the orthorombic cell
        side = [float(L) for L in fh.read().split()]

    with TrajectoryXYZ(file_inp) as th:
        th.add_callback(fix_missing_cell, side)

As a more permanent solution, you can define your own custom trajectory by subclassing ``TrajectoryXYZ``. First, parse the cell information during the initialization stage (``read_init()``).

.. code:: python

    from atooms.system import Cell
    from atooms.trajectory import TrajectoryXYZ

    class TrajectoryCustomXYZ(TrajectoryXYZ):

        def read_init(self):
            super().read_init()
            with open(self.filename + '.cell') as fh:
                self._side = [float(L) for L in fh.read().split()]

Then modify the ``read_sample()`` method, which reads a given frame of the trajectory.

.. code:: python

    def read_sample(self, frame):
        system = super().read_sample()
        system.cell = Cell(self._side)
        return system

Here we have assumed that the cell side is the same for all frames. The code would have to be adjusted to the more general case of a fluctuating cell.
