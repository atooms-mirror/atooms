import os
import shutil
import logging

# Logging facilities

LOGGER_NAME = 'atooms'
DEFAULT_LOGGING_FORMAT = '[%(levelname)s/%(processName)s] %(message)s'

_logger = None

def log_to_stderr(level=None):
    '''
    Turn on logging and add a handler which prints to stderr
    '''
    import logging

    logger = logging.get_logger(LOGGER_NAME)
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
    logging.info('found mpi4py %d %d' % (rank, size))
except:
    comm = None
    rank = 0
    size = 1
    logging.info('mpi4py not found')

def barrier():
    if size > 1:
        comm.barrier()

# Utility functions to mimic bash directory / file handling

def mkdir(d):
    if isinstance(d, str):
        dirs = [d]
    else:
        dirs = d

    for dd in dirs:
        try:
            os.makedirs(dd)
        except:
            pass

def rmd(files):
    try:
        shutil.rmtree(files)
    except:
        pass

def rmf(files):
    try:
        os.remove(files)
    except:
        pass


# Timer class, inspired by John Paulett's stopwatch

import time

class Timer(object):

    def __init__(self):
        self.__start_cpu = None
        self.__start_wall = None
        self.cpu_time = 0.0
        self.wall_time = 0.0
 
    def start(self):
        self.__start_cpu = self.__now_cpu()
        self.__start_wall = self.__now_wall()

    def stop(self):
        if self.__start_cpu is None:
            raise ValueError("Timer not started")
        self.cpu_time += self.__now_cpu() - self.__start_cpu
        self.wall_time += self.__now_wall() - self.__start_wall
   
    def __now_cpu(self):
        return time.time()

    def __now_wall(self):
        try:
            return MPI.Wtime()
        except:
            return 0.0
   

def clockit(func):
    """Function decorator that times the evaluation of *func* and prints the
    execution time.
    """
    def new(*args, **kw):
        t = Timer()
        t.start()
        retval = func(*args, **kw)
        t.stop()
        print '%s in %s' % (func.__name__, t)
        del t
        return retval
    return new


def fractional_slice(first, last, skip, n):
    """Return a slice assuming first or last are fractions of n, the length of the iterable,
    if first or last are in (0,1)"""
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
