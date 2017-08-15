#!/usr/bin/env python
#
# File: filter.py
#
# Created by: au447708 on 8/15/17
#
# Description:

import os
import cPickle
from config import cfg

def main():
    wdir = os.path.join(os.environ['WORK'], 'cdp')
    n = int(os.environ['SLURM_JOB_NUM_NODES'])
    p = int(os.environ['SLURM_PROCID'])
    t = int(os.environ['SLURM_NTASKS_PER_NODE'])
    steps = range(p, cfg['max_step']+1, n*t)
    for step in steps:
        try:
            ndir = os.path.join(wdir, 'step{}'.format(step),
                                'n{}'.format(cfg['min_for_cluster']))
            os.mkdir(ndir)
        except OSError:
            print('{} exists. Checking rotation file.'.format(ndir))
            if os.path.exists(ndir + '/rotations.pkl'):
                print('Skipping step {}.'.format(step))
                continue
        inrot = os.path.join(wdir, 'step{}'.format(step), 'rotations.pkl')
        try:
            fh = open(inrot)
        except IOError:
            print('Could not open {}. Skipping step {}'.format(inrot, step))
            continue
        else:
            with fh:
                pr = cPickle.load(fh)
        outrot = {p: pr[p] for p in pr if len(pr[p]) >= cfg['min_for_cluster']}
        outf = os.path.join(ndir, 'rotations.pkl')
        with open(outf, 'wb') as o:
            cPickle.dump(outrot, o, -1)


if __name__ == '__main__':
    main()
