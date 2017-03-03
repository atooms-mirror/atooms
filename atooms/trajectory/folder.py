# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Folder based trajectory.

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


class TrajectoryFolder(TrajectoryBase):

    def __init__(self, filename, mode='r', file_pattern='*', step_pattern='(\d*)'):
        TrajectoryBase.__init__(self, filename.rstrip('/'), mode)
        self.file_pattern = file_pattern
        self.step_pattern = step_pattern
        if mode == 'r':
            # See if trajectory is packed as a compressed tar file.
            # If so, configurations will be extracted inplace and deleted at the end.
            try:
                dirname = tempfile.mkdtemp()
                with tarfile.open(filename) as th:
                    th.extractall(path=dirname)
                    files = [os.path.join(dirname, f.name) for f in th.getmembers()]
                self.dirname = dirname
                self.archive = True
            except:
                if not os.path.isdir(filename):
                    raise IOError("Directory expected (%s)" % filename)
                self.dirname = filename
                self.archive = False
                files = glob.glob(os.path.join(self.dirname, self.file_pattern))

            self.files, self.steps = self._get_file_steps(files)

    def _get_step(self, fileinp):
        # Make sure we only test the basename (avoid metching
        # patterns in directory path)
        s = re.search(self.step_pattern, os.path.basename(fileinp))
        if s:
            step = int(s.group(1))
            return step
        else:
            raise ValueError('Could not find step')          

    def _get_file_steps(self, files):
        """Return a list of tuples (file, step). The step is extracted from
        the file path using `regexp`, which must contain one group for
        the step.
        """
        file_steps = []
        for i, f in enumerate(files):
            if os.path.isdir(f):
                continue
            try:
                step = self._get_step(f)
            except ValueError:
                step = i+1
            file_steps.append((f, step))
        file_steps.sort(key = lambda a : a[1])
        return [a[0] for a in file_steps], [a[1] for a in file_steps]

    def close(self):
        if self.archive:
            rmd(self.dirname)




