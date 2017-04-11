# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""
Folder based trajectory.

This is a base class for trajectories that stores configurations as
individual files in a directory. The `dirname` property stores the
directory path; an additional `files` property stores the list of
files sorted by step index.

It supports compressed archives. In this case the input filename is
the path to the archive, which will be decompressed in a temporary
folder.
"""

import os
import glob
import re
import tarfile
import tempfile

from atooms.utils import rmd
from atooms.trajectory.base import TrajectoryBase

# Helper functions

def init_folder(filename, file_pattern='*', step_pattern='(\d*)'):
    path = filename.rstrip('/')
    # See if trajectory is packed as a compressed tar file.
    # If so, configurations will be extracted inplace and deleted at the end.
    try:
        dirname = tempfile.mkdtemp()
        with tarfile.open(filename) as th:
            th.extractall(path=dirname)
            files = [os.path.join(dirname, f.name) for f in th.getmembers()]
        dirname = dirname
        archive = True
    except:
        if not os.path.isdir(filename):
            raise IOError("Directory expected (%s)" % filename)
        dirname = filename
        archive = False
        files = glob.glob(os.path.join(dirname, file_pattern))

    files, steps = _get_file_steps(files, step_pattern)
    return dirname, archive, files, steps

def _get_step(fileinp, step_pattern):
    # Make sure we only test the basename (avoid metching
    # patterns in directory path)
    s = re.search(step_pattern, os.path.basename(fileinp))
    if s:
        step = int(s.group(1))
        return step
    else:
        raise ValueError('Could not find step')          

def _get_file_steps(files, step_pattern):
    """Return a list of tuples (file, step). The step is extracted from
    the file path using `regexp`, which must contain one group for
    the step.
    """
    file_steps = []
    for i, f in enumerate(files):
        if os.path.isdir(f):
            continue
        try:
            step = _get_step(f, step_pattern)
        except ValueError:
            step = i+1
        file_steps.append((f, step))
    file_steps.sort(key = lambda a : a[1])
    return [a[0] for a in file_steps], [a[1] for a in file_steps]


# Classes

class TrajectoryFolder(TrajectoryBase):

    """Folder based trajectory."""

    def __init__(self, filename, mode='r', file_pattern='*', step_pattern='(\d*)'):
        TrajectoryBase.__init__(self, filename.rstrip('/'), mode)
        self.dirname, self.archive, self.files, self.steps = init_folder(filename, file_pattern, step_pattern)

    def close(self):
        if self.archive:
            rmd(self.dirname)


class Foldered(TrajectoryFolder):

    """Transform a file-based trajectory into folder-based one. Read-only."""

    def __init__(self, filename, mode='r', cls=None, file_pattern='*', step_pattern='(\d*)'):
        if mode != 'r':
            raise ValueError('Not ready for write mode')
        TrajectoryFolder.__init__(self, filename, mode)
        self._cls = cls
        if self.mode == 'r':
            self.dirname, self.archive, self.files, self.steps = init_folder(self.filename)

    def read_sample(self, sample):
        from atooms.trajectory import Trajectory
        with Trajectory(self.files[sample], fmt=self._cls) as th:
            return th.read_sample(0)

    def close(self):
        if self.archive:
            rmd(self.dirname)

        





