# This file is part of atooms
# Copyright 2010-2014, Daniele Coslovich

"""Folder based trajectory.

This is a base class for trajectories that stores configurations as
individual files in a directory. The `dirname` property stores the
directory path; an additional `files` property stores the list of
files sorted by step index.
"""

import os
import glob
import re

from atooms.trajectory.base import TrajectoryBase


class TrajectoryFolder(TrajectoryBase):

    def __init__(self, dirname, mode='r', file_pattern='*', step_pattern='(\d*)'):
        if not os.path.isdir(dirname):
            raise IOError("Directory expected (%s)" % dirname)
        TrajectoryBase.__init__(self, dirname, mode)
        self.dirname = dirname
        all_files = glob.glob(os.path.join(dirname, file_pattern))
        file_steps = self._get_file_steps(all_files, step_pattern)
        self.files = [a[0] for a in file_steps]
        self.steps = [a[1] for a in file_steps]

    def _get_file_steps(self, files, regexp):
        """Return a list of tuples (file, step). The step is extracted from
        the file path using `regexp`, which must contain one group for
        the step.
        """
        file_steps = []
        for i, f in enumerate(files):
            # Make sure we only test the basename (avoid metching
            # patterns in directory path)
            fbase = os.path.basename(f)
            s = re.search(regexp, fbase)
            if s:
                try:
                    step = int(s.group(1))
                except:
                    step = i+1
                file_steps.append((f, step))
        file_steps.sort(key = lambda a : a[1])
        return file_steps
