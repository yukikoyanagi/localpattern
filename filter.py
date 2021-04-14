#!/usr/bin/env python
#
# File: filter.py
#
# Description: Filter pattern data where the number of occurrences is less
# than the minimum specified in config.yaml.


import os.path
import argparse
import cPickle
import glob
from config import cfg

def get_rot(rot_txt):
    uv, phi = rot_txt.split(':')
    phi = float(phi)
    x, y, z = map(float, uv.split(','))
    return (x, y, z, phi)
    
            
def main(ind, outd, step):
    pr = {}

    for fn in glob.glob('{}/*.patrot'.format(ind)):
        with open(fn) as fh:
            for line in fh:
                s, pat, rots = line.split()
                s = int(s)
                if s==step:
                    rot = get_rot(rots)
                    try:
                        pr[pat].append(rot)
                    except KeyError:
                        pr[pat] = [rot]

    results = {p: pr[p] for p in pr if len(pr[p]) >= cfg['min_for_cluster']}
    outf = os.path.join(outd, 'step{}_rotations.pkl'.format(step))
    with open(outf, 'wb') as fh:
        cPickle.dump(results, fh, -1)
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_dir',
                        help='Directory with pattern-rotation files.')
    parser.add_argument('out_dir',
                        help='Output directory.')
    parser.add_argument('step', type=int,
                        help='Step number to process.')
    args = parser.parse_args()
    main(args.in_dir, args.out_dir, args.step)
