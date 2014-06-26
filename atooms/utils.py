import os
import shutil

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
