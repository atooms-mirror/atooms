# This file is part of atooms
# Copyright 2010-2017, Daniele Coslovich

"""PDB format, write-only."""

from .base import TrajectoryBase
from atooms.system.particle import Particle
from atooms.system.cell import Cell
from atooms.system import System


class TrajectoryPDB(TrajectoryBase):

    """PDB format (https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format))"""

    suffix = 'pdb'

    def __init__(self, filename, mode='r'):
        super(TrajectoryPDB, self).__init__(filename, mode)
        self._file = open(self.filename, self.mode)
        if mode == 'r':
            self._setup_index()

    def read_steps(self):
        pass

    def _setup_index(self):
        self._index_frame = []
        self.steps = []
        fh = self._file
        self._index_header_lines = 0
        header = True
        while True:
            line = fh.tell()
            data = fh.readline()
            # We break if file is over or we found an empty line
            if not data:
                break

            # Count number of header lines
            if header:
                self._index_header_lines += 1

            # Store lines of frame beginnings
            if data.startswith('MODEL'):
                # End of the header
                if header:
                    header = False
                self.steps.append(int(data.split()[-1]))
                self._index_frame.append(line)

    def write_system(self, system, step):
        cfg = ''
        cfg += 'MODEL%9i\n' % step
        fmt = 'CRYST1' + 3*'{:9.3f}' + 3*'{:7.2f}' + ' P 1           1\n'
        cfg += fmt.format(*(list(system.cell.side) + [90, 90, 90]))
        for p in system.particle:
            # If particle has a field property we dump it in the pdb file
            if hasattr(p, 'field'):
                x = float(p.field)
            else:
                x = 1.0
            data = [' '] * 80
            data[0:6] = 'HETATM'
            data[30:38] = '{:8.3f}'.format(p.position[0])
            data[38:46] = '{:8.3f}'.format(p.position[1])
            data[46:54] = '{:8.3f}'.format(p.position[2])
            data[60:66] = str(x)
            data[76:77] = p.species
            cfg += ''.join(data) + '\n'
        cfg += 'ENDMDL\n'
        self._file.write(cfg)

    def read_init(self):
        self._cell = None
        self._file.seek(0)
        for _ in range(self._index_header_lines):
            data = self._file.readline()
            if data.startswith('CRYST1'):
                Lx = float(data[6:15])
                Ly = float(data[15:24])
                Lz = float(data[24:33])
                self._cell = Cell([Lx, Ly, Lz])

    def read_system(self, frame):
        particle = []
        cell = self._cell

        # Read system at frame
        self._file.seek(self._index_frame[frame])
        _ = self._file.readline()
        while True:
            data = self._file.readline()
            if not data or data.startswith('ENDMDL'):
                break
            if data.startswith('CRYST1'):
                Lx = float(data[6:15])
                Ly = float(data[15:24])
                Lz = float(data[24:33])
                cell = Cell([Lx, Ly, Lz])
            elif data.startswith('HETATM') or data.startswith('ATOM'):
                atom = data[0:6]
                serial = data[6:11]
                name = data[12:16]
                altLoc = data[16]
                resName = data[17:20]
                chainID = data[21]
                resSeq = data[22:26]
                iCode = data[26]
                x = float(data[30:38])
                y = float(data[38:46])
                z = float(data[46:54])
                occupancy = data[54:60]
                tempFactor = data[60:66]
                try:
                    element = data[76:77]
                    charge = data[78:80]
                except IndexError:
                    element = name
                particle.append(Particle(position=[x, y, z], species=element))
        system = System()
        system.particle = particle
        system.cell = cell
        return system

    def close(self):
        self._file.close()
