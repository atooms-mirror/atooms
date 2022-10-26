import numpy
from .base import TrajectoryBase
from atooms.system import System, Particle


class TrajectoryCSV(TrajectoryBase):

    def __init__(self, filename, mode='r'):
        super(TrajectoryCSV, self).__init__(filename, mode)
        self._file = open(self.filename, self.mode)
        self.precision = 4
        self._ndim = 3
        self._sep = None
        self._order = 'C'
        if mode == 'r':
            self._setup()

    def _setup(self):
        self._file_index = []
        self._file.seek(0)
        while True:
            line = self._file.tell()
            data = self._file.readline()
            if not data:
                break
            else:
                self._file_index.append(line)
        self._file.seek(0)
        self.steps = range(len(self._file_index))

    def read_system(self, frame):
        self._file.seek(self._file_index[frame])
        data = self._file.readline().strip()
        pos = numpy.array([float(_) for _ in data.split(self._sep)])
        npart = len(pos) // self._ndim
        pos = pos.reshape((npart, self._ndim))
        system = System()
        for i in range(npart):
            #p = Particle(position=pos[i*self._ndim: (i+1)*self._ndim])
            p = Particle(position=pos[i, :])
            p.position_unfolded = p.position
            system.particle.append(p)
        return system

    def write_system(self, system, step):
        sep = ' '
        fmt = '{:.' + str(self.precision) + 'g}'
        numpy.set_string_function(lambda x: sep.join([fmt.format(xi) for xi in x]), False)
        self._file.write(str(system.dump('position', flat=True, order=self._order)) + '\n')
        numpy.set_string_function(None, False)

    def close(self):
        self._file.close()
