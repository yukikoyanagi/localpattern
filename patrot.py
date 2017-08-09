#!/usr/bin/env python
#
# File: patrot.py
#
# Description: Produce dict {pattern: [rotation]} for clustering. Output is
# written to a python pickle file.
#
# Author: Yuki Koyanagi
#
import os
import cPickle
from glob import glob
from collections import Counter

from config import cfg
import Protein
import Option


def getsteps(p):
    """Compute a list of steps to process on this core."""
    if p['abacus']:
        # n = int(os.environ['SLURM_JOB_NUM_NODES'])
        # return range(int(os.environ['SLURM_PROCID']),
        #              p['max_step'] + 1,
        #              p['cpu_per_node'] * n)
        fs = glob('/work/austmathjea/cdp/step*/rotations.pkl')
        e = [int(os.path.basename(os.path.dirname(f)).lstrip('step'))
             for f in fs]
        s = [i for i in range(792) if i not in e]
        try:
            l = [s[int(os.environ['SLURM_PROCID'])]]
        except IndexError:
            l = []
        return l
    else:
        return range(p['max_step'] + 1)


def filterpat(pats, minn):
    """
    Removes elements from the list which occur less than the given threshold.
    :param pats: List of elements to filter
    :param minn: Min. threshold
    :return: List of (unique) elements, all with occurrence >= minn
    """
    cnt = Counter(pats)

    return [k for k in cnt.keys() if cnt[k] >= minn]


def main():
    fs = glob(cfg['protdir'] + '/*.txt')
    steps = getsteps(cfg)

    for step in steps:
        outfile = os.path.join(cfg['outdir'], 'step{}'.format(step),
                               'rotations.pkl')
        if os.path.exists(outfile):
            print('{} already exists.'.format(outfile))
            continue
        of = os.path.join(cfg['optsdir'], 'step{}_opts'.format(step))
        opt = Option.Option(of)
        patrot = {}
        for f in fs:
            protid = os.path.splitext(os.path.basename(f))[0]
            tertf = os.path.join(cfg['tertdir'], protid + '.txt')
            if not os.path.exists(tertf):
                # Corresponding tertiary file may not be present...
                tertf = None
            prot = Protein.Protein(protid)
            prot.fromfiles(f, tertf)
            for bond in prot.hbonds:
                try:
                    patrot[str(prot.findpattern2(bond, opt))].append(
                        bond.rotation)
                except KeyError:
                    patrot[str(prot.findpattern2(bond, opt))] = [bond.rotation]

        with open(outfile, 'wb') as o:
            cPickle.dump(patrot, o)


if __name__ == '__main__':
    main()
