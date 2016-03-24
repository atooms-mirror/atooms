#!/usr/bin/env python

import os
import sys
import h5py
import numpy
#from pyatooms.trajectory import trajectory
from pyutils.utils import mkdir

def pratio(vec):
    pra = numpy.sum([numpy.dot(ri, ri)**2 for ri in vec])
    return 1.0 / (pra * vec.shape[0])

def main(f):
    fh = h5py.File(f)
    samples = [d[:] for d in fh["trajectory/realtime/sampleindex"]]
    steps = [d[0] for d in fh["trajectory/realtime/stepindex"].values()]
    try:
        modes = [fh["trajectory/normalmodes/eigenvectors/index/%s" % s][:] for s in samples]
    except:
        raise IOError('file does not contain eigenvectors')

    for s, step, mm in zip(samples[:], steps[:], modes[:]):
        l = "trajectory/normalmodes/eigenvalues/%s" % s
        out_dir = f + '.dump/step_%s' % step
        mkdir(out_dir)
        file_val = '%s/eigenvalues.txt' % out_dir
        fh_val = open(file_val, 'w')
        fh_val.write('# file eigenvalue omega pratio\n')
        for mi in mm:
            g = "trajectory/normalmodes/eigenvectors/vector/%s_mode_%05d" % (s, mi)
            file_eigvec = '%s/eigenvector_%d.txt' % (out_dir, mi)
            lam, pr, vec = fh[l][mi], pratio(fh[g][:]), fh[g][:]
            fh_val.write('%s %g %g %g\n' % (os.path.basename(file_eigvec), lam, lam**0.5, pr))
            with open(file_eigvec, 'w') as fh_vec:
                fh_vec.write('%d\n' % len(vec))
                fh_vec.write('w:%f pratio:%f\n' % (lam**0.5, pr))
                for x in vec:
                    fh_vec.write('%g %g %g\n' % (x[0], x[1], x[2]))
        fh_val.close()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='verbose', action='store_true', help='verbose output')
    parser.add_argument(nargs='+', dest='files',type=str, help='input files')
    args = parser.parse_args()

    for f in args.files:
        main(f)
    
