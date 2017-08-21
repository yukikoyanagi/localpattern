#!/usr/bin/env python
#
# File: patrot_prll.py
#
# Created by: au447708 on 8/10/17
#
# Description: Produce dict {pattern: [rotation]} for clustering, using
# ParallelPython. Output is written to a python pickle file.

import os
import cPickle
import uuid
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


def findpatterns(step, data):
    """
    Find local patterns for the given step & data.
    :param step: int; step number
    :param data: list of prot file names to process.
    :return: dict of {pattern: list of rotations}
    """
    import os
    import uuid
    import cPickle
    import Protein
    import Option
    from config import cfg

    sdir = os.environ['SCRATCH']
    of = os.path.join(cfg['optsdir'], 'step{}_opts'.format(step))
    opt = Option.Option(of)
    patrot = {}
    for f in data:
        protid = os.path.splitext(os.path.basename(f))[0]
        tertf = os.path.join(cfg['tertdir'], protid + '.txt')
        if (not os.path.exists(tertf)) or (cfg['max_tbond_level'] == -1):
            # Corresponding tertiary file may not be present...
            tertf = None
        prot = Protein.Protein(protid)
        prot.fromfiles(f, tertf)
        for bond in prot.hbonds:
            pat = str(prot.findpattern2(bond, opt))
            if pat in patrot:
                patrot[pat].append(tuple(bond.rotation))
            else:
                patrot[pat] = [tuple(bond.rotation)]
        del prot

    outfile = os.path.join(sdir,
                           'step{}_{}.pkl'.format(step, uuid.uuid4().hex))
    with open(outfile, 'wb') as o:
        cPickle.dump(patrot, o, -1)


def collect(step):
    """
    Collect intermediate rotation files to a single file per step.
    :param step: int; step number
    :return:
    """
    import os
    import cPickle
    from glob import glob
    from config import cfg

    results = {}
    tfs = glob(os.environ['SCRATCH'] + '/step{}_*.pkl'.format(step))

    for tf in tfs:
        with open(tf, 'rb') as fh:
            d = cPickle.load(fh)
        for k in d:
            if k in results:
                results[k] += d[k]
            else:
                results[k] = d[k]

    outfile = os.path.join(cfg['outdir'], 'step{}'.format(step),
                           'rotations.pkl')
    with open(outfile, 'wb') as o:
        cPickle.dump(results, o, -1)

    # Clean up temp files
    for tf in tfs:
        os.remove(tf)


def main():
    steps = range(cfg['max_step'] + 1)
    fs = glob(cfg['protdir'] + '/*.txt')
    subsets = []
    for k, g in groupby(enumerate(fs), key=lambda x: x[0] / 2000):
        subsets.append([f[1] for f in g])

    submitted = []

    for step in steps:
        outfile = os.path.join(cfg['outdir'], 'step{}'.format(step),
                               'rotations.pkl')
        if os.path.exists(outfile):
            print('{} already exists.'.format(outfile))
            continue

        for subset in subsets:
            job_server.submit(findpatterns,
                              args=(step, subset),
                              group='step{}'.format(step))
        submitted.append(step)

    job_server.wait()

    for step in submitted:
        job_server.submit(collect, args=(step,))

    job_server.wait()

if __name__ == '__main__':
    main()