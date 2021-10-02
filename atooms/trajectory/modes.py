import os
import numpy
from numpy import copysign
from atooms.system import System
from atooms.trajectory.hdf5 import TrajectoryHDF5


class TrajectoryModes(TrajectoryHDF5):

    def read_system(self, frame):
        system = System()
        s = self._file["trajectory/realtime/sampleindex"].keys()[frame]
        system.eigenvalues = self._file["trajectory/normalmodes/eigenvalues/%s" % s][:]
        system.eigenfreq = numpy.array([copysign(abs(x)**0.5, x) for x in system.eigenvalues])
        mode_idx = self._file["trajectory/normalmodes/eigenvectors/index/%s" % s][:]
        eigv = {}
        for idx in mode_idx:
            # Modes are F-indexed
            omega = system.eigenfreq[idx-1]
            vect = self._file["trajectory/normalmodes/eigenvectors/vector/%s_mode_%05d" % (s, idx)][:]
            eigv[omega] = vect
        system.eigenvectors = eigv

        # Add positions
        parent = os.path.splitext(self._file.filename)[0]
        step = self.steps[frame]
        with TrajectoryHDF5(parent) as th:
            # The step should be there
            parent_frame = th.steps.index(step)
            system.particle = th[parent_frame].particle
            system.cell = th[parent_frame].cell
        return system

# def test(f):
#     with TrajectoryModes('/tmp/config.h5.sad.inm') as inm:
#         print inm[0]
#         print inm[0].particle
#         print inm[0].eigenvalues[:4] # array, full
#         print inm[0].eigenvectors.values()[:4] # dict, what if we have two identical eigenvalues? or list of arrays but we need the list of w

# test('')


if __name__ == '__main__':
    import sys
    with TrajectoryModes(sys.argv[1]) as inm:
        for i in range(len(inm.steps)):
            w = inm[i].eigenfreq[1]
            # dict, what if we have two identical eigenvalues? or list of arrays but we need the list of w
            print(w, inm[i].particle[0].position, inm[i].eigenvectors[w][0])
