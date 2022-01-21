import os
import bz2
import numpy
from xml.etree import ElementTree
from atooms.trajectory.base import TrajectoryBase, SuperTrajectory
from atooms.system import Particle, Cell, System


class _TrajectoryDynamO(TrajectoryBase):

    """DynamO trajectory format (https://www.dynamomd.com/index.php/tutorial3)"""

    suffix = 'xml'

    def __init__(self, filename, mode='r'):
        super(_TrajectoryDynamO, self).__init__(filename, mode)
        assert mode == 'r', 'cannot write yet'
        self._filename = filename
        if (os.path.splitext(filename)[1][1:].strip() == "bz2"):
            self._file = bz2.BZ2File(filename)
        else:
            self._file = open(filename)
        self.variables = ['particle.species', 'particle.mass',
                          'particle.position', 'particle.velocity']

    def read_timestep(self):
        return 1e-10

    def read_steps(self):
        import re
        base = os.path.basename(self._filename)
        # If the file is named with Snapshot.[0-9]*e.xml* naming
        # scheme then we extract the time from the corresponding
        # output file and convert to a step assuming that frames are
        # stored at regular time intervals. Otherwise the step is 0.
        if base.startswith('Snapshot'):
            # Get the current frame time from output file
            output_file = re.sub('([0-9].)', 'output.\\1', self._filename)
            if (os.path.splitext(output_file)[1][1:].strip() == "bz2"):
                fh = bz2.BZ2File(output_file)
            else:
                fh = open(output_file)
            tree = ElementTree.parse(fh)
            root = tree.getroot()
            time = root.find('Misc').find('Duration').attrib['Time']
            step = round(float(time) / self.timestep)
            fh.close()
        else:
            step = 0
        return [step]

    def read_system(self, frame):
        assert frame == 0
        tree = ElementTree.parse(self._file)
        root = tree.getroot()

        # Cell
        _cell = root.find("Simulation").find('SimulationSize')
        side = numpy.array([_cell.attrib['x'], _cell.attrib['y'], _cell.attrib['z']])

        # Coordinates
        position, velocity = [], []
        for _particle in root.findall('.//Pt'):
            pos = _particle.find('P')
            vel = _particle.find('V')
            position.append([pos.attrib['x'], pos.attrib['y'], pos.attrib['z']])
            velocity.append([vel.attrib['x'], vel.attrib['y'], vel.attrib['z']])

        # Species
        _species = root.find("Simulation").find("Genus").findall('Species')
        if len(_species) == 1:
            species, mass = len(position) * ['A'], len(position) * [1.0]
        else:
            species, mass = [], []
            for sp in _species:
                assert sp.find('IDRange').attrib['Type'] == 'Ranged'
                begin = int(sp.find('IDRange').attrib['Start'])
                end = int(sp.find('IDRange').attrib['End'])
                for _ in range(begin, end+1):
                    species.append(sp.attrib['Name'])
                    mass.append(sp.attrib['Mass'])

        # Pack everything into a System
        cell = Cell(side)
        particle = [Particle(species=s, mass=float(m),
                             position=numpy.array(p, dtype=float),
                             velocity=numpy.array(v, dtype=float))
                    for p, m, s, v in zip(position, mass, species, velocity)]
        return System(particle, cell)

    def close(self):
        self._file.close()


class TrajectoryDynamO(TrajectoryBase):

    """DynamO trajectory format (https://www.dynamomd.com/index.php/tutorial3)"""

    def __init__(self, path, mode='r'):
        import os
        import glob

        super(TrajectoryDynamO, self).__init__(path, mode)
        if os.path.isdir(path):
            files = glob.glob(path + '/Snapshot.[0-9]*.*')
            self._th = SuperTrajectory(files, _TrajectoryDynamO)
        else:
            self._th = _TrajectoryDynamO(path)
        self.steps = self._th.steps

    def read_timestep(self):
        return self._th.read_timestep()

    def read_system(self, frame):
        return self._th.read_system(frame)
