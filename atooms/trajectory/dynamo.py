"""DynamO trajectory format"""

import os
import bz2
import numpy
from xml.etree import ElementTree
from atooms.trajectory.base import TrajectoryBase
from atooms.system import Particle, Cell, System


class TrajectoryDynamO(TrajectoryBase):

    suffix = 'xml'

    def __init__(self, filename, mode='r'):
        super(TrajectoryDynamO, self).__init__(filename, mode)
        assert mode == 'r', 'cannot write yet'
        self._filename = filename
        if (os.path.splitext(filename)[1][1:].strip() == "bz2"):
            self._file = bz2.BZ2File(filename)
        else:
            self._file = open(filename)
        self.variables = ['particle.species', 'particle.mass',
                          'particle.position', 'particle.velocity']

    def read_steps(self):
        # Quick and dirty (assume Snapshot.[0-9]*e.xml* naming scheme)
        # TODO: look up output file if it exists
        step = 0
        base = os.path.basename(self._filename)
        if base.startswith('Snapshot'):
            step = int(base.split('.')[1][:-1])
        return [step]

    def read_system(self, frame):
        assert frame == 0
        tree = ElementTree.parse(self._file)
        root = tree.getroot()

        # Cell
        _cell = root.find("Simulation").find('SimulationSize')
        side = numpy.array([_cell.attrib['x'], _cell.attrib['y'], _cell.attrib['z']])

        # Species
        species, mass = [], []
        _species = root.find("Simulation").find("Genus").findall('Species')
        for sp in _species:
            assert sp.find('IDRange').attrib['Type'] == 'Ranged'
            begin = int(sp.find('IDRange').attrib['Start'])
            end = int(sp.find('IDRange').attrib['End'])
            for _ in range(begin, end+1):
                species.append(sp.attrib['Name'])
                mass.append(sp.attrib['Mass'])

        # Coordinates
        position, velocity = [], []
        for _particle in root.findall('.//Pt'):
            pos = _particle.find('P')
            vel = _particle.find('V')
            position.append([pos.attrib['x'], pos.attrib['y'], pos.attrib['z']])
            velocity.append([vel.attrib['x'], vel.attrib['y'], vel.attrib['z']])

        cell = Cell(side)
        particle = [Particle(species=s, mass=float(m),
                             position=numpy.array(p, dtype=float),
                             velocity=numpy.array(v, dtype=float))
                    for p, m, s, v in zip(position, mass, species, velocity)]
        return System(particle, cell)

    def close(self):
        self._file.close()
