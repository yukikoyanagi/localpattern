#!/usr/bin/env python
#
# File: clstr2.py
#
# Description: Run cluster analysis on the rotation data.
#

import os, cPickle, subprocess, argparse, os.path, uuid

from config import cfg

minn = cfg['min_for_cluster']
maxstep = cfg['max_step']


def runclstr(pattern, data, step, outd, dbg):
    rscript = 'Rasmus2.r'
    imprv = 'improve_modebox.pl'
    rotf = os.path.join(outd, uuid.uuid4().hex)

    with open(rotf, 'w') as f:
        f.write('\n'.join(['\t'.join(map(str, d)) for d in data]))

    try:
        subprocess.check_output(['Rscript', rscript, rotf],
                                stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        print 'Error while processing step{}/{}'.format(step, pattern)
        print cpe.output
        return None

    # sumf = rotf + '_Summary.txt'

    cmd = ' '.join(['perl', imprv, rotf])
    os.system(cmd)

    sumf = rotf + '_Summary2.txt'

    with open(sumf) as f:
        res = f.read()

    if not dbg:
        os.remove(rotf)
        bf = rotf + '_BoxId.txt'
        os.remove(bf)
        sf = rotf + '_Summary.txt'
        os.remove(sf)
        os.remove(sumf)

    return res



def main(inf, dbg):
    outd = os.path.dirname(inf)
    step = int(os.path.basename(inf).split('_')[0][4:])
    d = {}
    with open(inf, 'rb') as fh:
        data = cPickle.load(fh)
    for pattern in data:
        # If R script has error runclstr returns None
        summ = runclstr(pattern, data[pattern], step, outd, dbg)
        if summ:
            d[pattern] = summ
    outf = os.path.join(outd, 'step{}_summary.pkl'.format(step))
    with open(outf, 'wb') as fh:
        cPickle.dump(d, fh)



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file',
                        help='A pickle file with pattern-rotation data, '
                        'named step##_rotations.pkl.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Keep all intermediate files.')
    args = parser.parse_args()
    main(args.in_file, args.debug)
