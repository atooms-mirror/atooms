import os
import shutil

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

