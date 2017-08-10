#!/usr/bin/env python
#
# File: patrot_prll.py
#
# Created by: au447708 on 8/10/17
#
# Description: Produce dict {pattern: [rotation]} for clustering, using
# ParallelPython. Output is written to a python pickle file.

import os
import cPickle as pickle
from glob import glob
from itertools import groupby

import pp

from config import cfg
import Protein
import Option

# Set up servers
ppservers = open('/tmp/nodelist').read().strip().split()
ppservers = tuple(pp + ':2048' for pp in ppservers)
job_server = pp.Server(0, ppservers=ppservers)


def findpatterns(step, cfg, data):
    """
    Find local patterns for the given step & data.
    :param step: int; step number
    :param cfg: config.cfg object
    :param data: list of prot file names to process.
    :return: dict of {pattern: list of rotations}
    """
    of = os.path.join(cfg['optsdir'], 'step{}_opts'.format(step))
    opt = Option.Option(of)
    patrot = {}
    for f in data:
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

    return patrot


def main():
    steps = range(cfg['max_step'] + 1)
    fs = glob(cfg['protdir'] + '/*.txt')
    subsets = []
    for k, g in groupby(enumerate(fs), key=lambda x: x[0] / 200):
        subsets.append([f[1] for f in g])

    jobs = []
    submitted = []

    for step in steps:
        outfile = os.path.join(cfg['outdir'], 'step{}'.format(step),
                               'rotations.pkl')
        if os.path.exists(outfile):
            print('{} already exists.'.format(outfile))
            continue

        for subset in subsets:
            job = job_server.submit(findpatterns,
                                    (step, cfg, subset),
                                    (),
                                    ("os", "Protein", "Option"),
                                    None,
                                    (),
                                    step)
            jobs.append(job)

        submitted.append(step)

    for step in submitted:
        job_server.wait(step)
        results = {}
        for job in [task for task in jobs if task.group == step]:
            d = job()
            results.update(d)

        outfile = os.path.join(cfg['outdir'], 'step{}'.format(step),
                               'rotations.pkl')
        with open(outfile, 'wb') as o:
            pickle.dump(results, o)


if __name__ == '__main__':
    main()