# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Linked cells to compute neighbors efficiently [*beta*]."""

from collections import defaultdict

def _pbc(t, N):
    # TODO: fix undefined NN
    for i in range(len(t)):
        if t[i] >= N[i]:
            t[i] -= NN
        elif t[i] < 0:
            t[i] += NN
    return t

def _pbc_all(l, N):
    return [map(_pbc, t) for t in l]

# TODO: define iterator over cells

class _LinkedCell(object):

    def __init__(self):
        self._particle_in_cell = None
        self.__all_cells = None
        self.__ghost_cells = None

    @property
    def _all_cells(self):
        if self.__all_cells is None:
            self.__all_cells = []
            for ix in range(self.n_cell[0]):
                for iy in range(self.n_cell[1]):
                    for iz in range(self.n_cell[2]):
                        self.__all_cells.append((ix, iy, iz))

        return self.__all_cells

    @property
    def _ghost_cells(self):
        if self.__ghost_cells is None:
            self.__ghost_cells = []
            for ix in range(self.n_cell[0]):
                for iy in range(self.n_cell[1]):
                    for iz in range(self.n_cell[2]):
                        if ix == -1 or ix == self.n_cell[0] or \
                           iy == -1 or iy == self.n_cell[1] or \
                           iz == -1 or iz == self.n_cell[2]:
                            self.__ghost_cells.append((ix, iy, iz))

        return self.__ghost_cells

    def _map_netwon(self):
        self._neigh_cell = {}
        for ix in range(self.n_cell[0]):
            for iy in range(self.n_cell[1]):
                for iz in range(self.n_cell[2]):
                    self._neigh_cell[(ix, iy, iz)] = \
                        [(ix+1, iy, iz  ), (ix+1, iy+1, iz  ), (ix  , iy+1, iz  ), (ix-1, iy+1, iz  ),
                         (ix+1, iy, iz-1), (ix+1, iy+1, iz-1), (ix  , iy+1, iz-1), (ix-1, iy+1, iz-1),
                         (ix+1, iy, iz+1), (ix+1, iy+1, iz+1), (ix  , iy+1, iz+1), (ix-1, iy+1, iz+1),
                         (ix  , iy, iz+1)]

        # Use minimum image convention and transform keys into numpy arrays
        for i in self._neigh_cell:
            for j in range(len(self._neigh_cell[i])):
                if self._neigh_cell[i][j] in self._ghost_cells:
                    self._neigh_cell[i][j] = _pbc(self._neigh_cell[i][j], self.n_cell)

    def adjust(self, npart, box, rcut):
        self.box = box
        self.hbox = box / 2
        self.n_cell = (box / max(rcut)).astype(int)
        self._map_netwon()

    def _index(self, pos, box):
        #return tuple(((pos + self.hbox) / box))
        return (pos + self.hbox) / box

    def compute(self, pos, box):
        # We only need positions here but how can we be sure that
        # this is the same set of particles we use when retrieving
        # the neighbours? We should keep a reference.
        self._particle_in_cell = {}
        # TODO: box should be cell
        index = self._index(pos, box)
        self._particle_in_cell = defaultdict(list)
        for ipart, i in enumerate(index):
            self._particle_in_cell[tuple(i)].append(ipart)

    def neighbours(self, pos):
        i = self._index(pos, self.box)
        neigh = []
        for j in self._neigh_cell[i]:
            neigh += self._particle_in_cell[j]
        return neigh
