"""Helper functions"""

import os
import sys
import shutil
import time


# Logging facilities

LOGGER_NAME = 'atooms'
DEFAULT_LOGGING_FORMAT = '[%(levelname)s/%(processName)s] %(message)s'

_logger = None

# We define the logging handler here to avoid "No handler found" warnings.
# Client classes should use this instead of logging.NullHandler
import logging
try:
    from logging import NullHandler
except ImportError:
    # Python <= 2.6
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass


def log_to_stderr(level=None):
    """
    Turn on logging and add a handler which prints to stderr
    """
    logger = logging.getLogger(LOGGER_NAME)
    formatter = logging.Formatter(DEFAULT_LOGGING_FORMAT)
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if level:
        logger.setLevel(level)
    return _logger


# Parallel environment

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    sys_excepthook = sys.excepthook

    def mpi_excepthook(v, t, tb):
        sys_excepthook(v, t, tb)
        MPI.COMM_WORLD.Abort(1)

    sys.excepthook = mpi_excepthook
except:
    comm = None
    rank = 0
    size = 1


def barrier():
    """Syncronize processes in parallel."""
    if size > 1:
        comm.barrier()


# Utility functions to mimic bash directory / file handling

def mkdir(dirname):
    """
    Create a directory `dirname` or a list `dirname` of directories,
    silently ignoring existing directories.

    This is just a wrapper to `os.makedirs`. All intermediate
    subdirectories are created as needed.
    """
    if dirname is None:
        return
    if isinstance(dirname, str):
        dirs = [dirname]
    else:
        dirs = dirname

    for dd in dirs:
        try:
            os.makedirs(dd)
        except OSError:
            pass


def rmd(files):
    """Totally silent wrapper to shutil.rmtree."""
    try:
        shutil.rmtree(files)
    except:
        pass


def rmf(files):
    """
    Remove `files` without complaining.

    The variable `files` can be a list or tuple of paths or a single
    string parseable by glob.glob().
    """
    import glob
    try:
        # This a single pattern
        for pathname in glob.glob(files):
            try:
                os.remove(pathname)
            except OSError:
                # File does not exists or it is a folder
                pass
    except (TypeError, AttributeError):
        # This is a list
        for pathname in files:
            try:
                os.remove(pathname)
            except OSError:
                # File does not exists or it is a folder
                pass


def cp(finp, fout):
    """
    Copy `finp` to `fout`.

    Wrapper to shutil.copy().
    """
    # Avoid overwriting file
    if finp == fout:
        return
    shutil.copy(finp, fout)

    
def wget(url, output_dir):
    """
    Python implementation of bash wget
    """
    import sys
    import os
    import shutil
    try:
        from urllib.request import urlopen # Python 3
    except ImportError:
        from urllib2 import urlopen # Python 2

    basename = os.path.basename(url)
    output_file = os.path.join(output_dir, basename)
    response = urlopen(url)
    length = 16*1024
    with open(output_file, 'wb') as fh:
        shutil.copyfileobj(response, fh, length)

# Alias of wget
download = wget


# Timings

class Timer(object):

    """Timer class inspired by John Paulett's stopwatch class."""

    def __init__(self):
        self.__start_cpu = None
        self.__start_wall = None
        self.cpu_time = 0.0
        self.wall_time = 0.0
        try:
            self._wall_time_func = MPI.Wtime
        except:
            self._wall_time_func = time.time

    def start(self):
        self.__start_cpu = self.__now_cpu()
        self.__start_wall = self.__now_wall()

    def stop(self):
        if self.__start_cpu is None:
            raise ValueError("Timer not started")
        self.cpu_time += self.__now_cpu() - self.__start_cpu
        self.wall_time += self.__now_wall() - self.__start_wall

    def __now_cpu(self):
        try:
            return time.clock()
        except AttributeError:
            return time.process_time()

    def __now_wall(self):
        return self._wall_time_func()


def clockit(func):
    """
    Function decorator that times the evaluation of `func` and prints
    the execution time.
    """
    def new(*args, **kw):
        t = Timer()
        t.start()
        retval = func(*args, **kw)
        t.stop()
        print('%s in %s' % (func.__name__, t))
        del t
        return retval
    return new


def fractional_slice(first, last, skip, n):
    """
    Return a slice assuming `first` or `last` are fractions of `n`,
    the length of the iterable, if `first` or `last` are in (0,1)
    """
    # We use an implicit convention here:
    # If first or last are in (0,1) then they are considered as fractions of the iterable
    # otherwise they are integer indexes. Note the explicit int() cast in the latter case.
    if first is not None:
        if first > 0 and first < 1:
            first = int(first * n)
        else:
            first = int(first)

    if last is not None:
        if last > 0 and last < 1:
            last = int(last * n)
        else:
            last = int(last)

    return slice(first, last, skip)


def add_first_last_skip(parser, what=None):
    """
    Add first, last, skip arguments to ArgumentParser object.

    Compatible with fractional_slice(). Convenience function for
    analysis scripts.
    """
    if what is None:
        what = ['first', 'last', 'skip']
    if 'first' in what:
        parser.add_argument('-f', '--first', dest='first', type=float, default=None, help='first cfg (accepts fractions)')
    if 'last' in what:
        parser.add_argument('-l', '--last', dest='last', type=float, default=None, help='last cfg (accepts fractions)')
    if 'skip' in what:
        parser.add_argument('-s', '--skip', dest='skip', type=int, default=1, help='interval between cfg')
    return parser


# Logging facilities

class _ParallelFilter(logging.Filter):
    def filter(self, rec):
        if hasattr(rec, 'rank'):
            if rec.rank == 'all':
                return True
            else:
                return rank == rec.rank
        else:
            return rank == 0


class _MyFormatter(logging.Formatter):
    def format(self, record):
        if record.levelname in ['WARNING', 'ERROR']:
            return '# ' + record.levelname + ' ' + record.msg % record.args
        else:
            return '# ' + record.msg % record.args


def setup_logging(name=None, level=40, filename=None, update=False):
    """Logging API."""
    if name is None:
        log = logging.getLogger()
    else:
        log = logging.getLogger(name)

    if update:
        # We only update the level of the logger
        log.setLevel(level)
        return
    else:
        # The logger should always pass messages to all handlers
        current_level = log.getEffectiveLevel()
        log.setLevel(min(level, current_level))

    formatter = _MyFormatter()
    if filename is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(filename)

    handler.setFormatter(formatter)
    # From the doc: "Note that filters attached to handlers are
    # consulted before an event is emitted by the handler, whereas
    # filters attached to loggers are consulted whenever an event is
    # logged (using debug(), info(), etc.), before sending an event to
    # handlers. This means that events which have been generated by
    # descendant loggers will not be filtered by a logger filter
    # setting, unless the filter has also been applied to those
    # descendant loggers."
    handler.addFilter(_ParallelFilter())
    handler.setLevel(level)
    if update:
        log.handlers(handler)
    else:
        log.addHandler(handler)

    return log


def tipify(s):
    """
    Convert a string into the best matching type.

    Example:
    -------
        2 -> int
        2.32 -> float
        text -> str

    The only risk is if a variable is required to be float,
    but is passed without dot.

    Tests:
    -----
        print type(tipify("2.0")) is float
        print type(tipify("2")) is int
        print type(tipify("t2")) is str
        print map(tipify, ["2.0", "2"])
    """
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


# Miscellaneous

def __header_dict(line):
    # Array entry have comma separated elements, split them into lists
    params = {}
    for key, value in [d.split('=') for d in line.split()]:
        params[key] = value
    return params


def report_parameters(params, fileout, version, comment=''):
    """Report parameters."""
    maxlen = max([len(key) for key in params])
    fmt = comment + '%-' + str(maxlen) + 's = %s\n'
    txt = ""
    txt += fmt % ('version', version)
    for key in sorted(params.keys()):
        txt += fmt % (key, params[key])
    if fileout is not None:
        with open(fileout, 'w') as fh:
            fh.write(txt)
    return txt

def report_command(cmd, params, main, fileout):
    """Report command line options."""
    txt = cmd + ' \\\n'
    for key in sorted(params.keys()):
        if key in main:
            continue
        flag = key.replace('_', '-')
        value = params[key]
        if value is None or value is False:
            continue
        if value is True:
            value = ''
        txt += ' --%s %s \\\n' % (flag, value)
    txt += ' '.join([params[key] for key in main])
    if fileout is not None:
        with open(fileout, 'w') as fh:
            fh.write(txt)
    return txt


class OrderedSet(object):

    """
    Simple class to store an ordered set of items.

    It does not try to reproduce either set or list interface. It just
    provides a simple interface to deal with the set of distinct
    chemical species in a system as it gets populated e.g. when
    reading a trajectory. This covers the use case of grand-canonical
    simulations.

    Example:
    -------

        particle = [Particle(species='A'), Particle(species='C')]
        periodic_table = OrderedSet()
        periodic_table.update([p.species for p in particle])
        particle = [Particle(species='A'), Particle(species='D')]
        periodic_table.update([p.species for p in particle])
        print periodic_table.index('C')
        print periodic_table[0]
        print periodic_table
    """

    def __init__(self):
        self.items = []

    def __repr__(self):
        return repr(self.items)

    def __iter__(self):
        return self.items.__iter__()

    def __getitem__(self, key):
        return self.items[key]

    def __setitem__(self, key, item):
        self.items.__setitem__
        return self.items[key]

    def update(self, items):
        sort_needed = False
        for item in items:
            if item not in self.items:
                self.items.append(item)
                sort_needed = True
        if sort_needed:
            self.items.sort()

    def index(self, item):
        try:
            return self.items.index(item)
        except ValueError:
            raise ValueError('item %s not in %s' % (item, self))

def is_array(arg):
    """
    Tell if `arg` is iterable but not string, morally an array.

    Source: https://stackoverflow.com/questions/1055360/how-to-tell-a-variable-is-iterable-but-not-a-string/44328500#44328500

    Examples:
    ---------
    # non-string iterables    
    assert iterable(("f", "f"))    # tuple
    assert iterable(["f", "f"])    # list
    assert iterable(iter("ff"))    # iterator
    assert iterable(range(44))     # generator
    assert iterable(b"ff")         # bytes (Python 2 calls this a string)

    # strings or non-iterables
    assert not iterable(u"ff")     # string
    assert not iterable(44)        # integer
    assert not iterable(iterable)  # function
    """
    import collections
    import six
    return isinstance(arg, collections.Iterable) and \
        not isinstance(arg, six.string_types)

def is_array(arg):
    """
    Non-python guess whether argument is an array.

    It will match list, tuple and numpy array but fail with iterators / generators.
    """
    import numpy
    return isinstance(arg, numpy.ndarray) or \
        isinstance(arg, list) or \
        isinstance(arg, tuple)
