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
from collections import defaultdict
from config import cfg

def get_rot(rot_txt):
    uv, phi = rot_txt.split(':')
    phi = float(phi)
    x, y, z = map(float, uv.split(','))
    return (x, y, z, phi)
    
            
def main(ind, outd, step, lf):
    pids_to_process = []
    if lf:
        with open(lf) as fh:
            pids_to_process = fh.read().splitlines()
    
    pr = defaultdict(list)

    for fn in glob.glob('{}/*.patrot'.format(ind)):
        if pids_to_process:
            pid = os.path.basename(fn).split('.')[0]
            if pid not in pids_to_process:
                continue
        with open(fn) as fh:
            for line in fh:
                s, pat, rots = line.split()
                s = int(s)
                if s==step:
                    rot = get_rot(rots)
                    pr[pat].append(rot)

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
    parser.add_argument('-l', '--list_file',
                        help='A file listing PIDs (one per line) to use. '
                        'Any file in in_dir that is not found in the list '
                        'is ignored.')
    args = parser.parse_args()
    main(args.in_dir, args.out_dir, args.step, args.list_file)
