#!/usr/bin/env python
#
# File: pred_prll.py
#
# Created by: au447708 on 8/14/17
#
# Description: Predict H-bond rotations. This is *almost* a clone of
# cdp.predict. Uses parallelpython (http://www.parallelpython.com).
# Requires: conv.py, Protein.py, Option.py, config.py.
#

import argparse
import os
import glob
import math
import cPickle
import re
import pp

import conv
from config import cfg


# Use the following pattern to construct a list of assess files
afpattern = '/work/austmathjea/cdp/step*/n{}/step*_assess'.format(
    cfg['min_for_cluster']
)

# Set up servers
ppservers = open('/tmp/nodelist').read().strip().split()
if len(ppservers) > 1:
    ppservers = tuple(pp + ':2048' for pp in ppservers)
    job_server = pp.Server(0, ppservers=ppservers, socket_timeout=10800)
else:
    job_server = pp.Server()


def runstep(step, pdir, af):
    '''
    Runs predictions for the given step. Returns a dict of
    {(protein, lineno): value}, value is a tuple containing
    pattern: string
    mode: 3-tuple of mode box id's
    ev: 4-tuple of rotation vector; angle, x, y, z
    score: score
    '''
    import os
    import shutil
    import cPickle
    import glob
    import Protein
    import Option
    from config import cfg

    # SLURM-specific environment var's
    try:
        os.environ['LOCALSCRATCH']
    except KeyError:
        os.environ['LOCALSCRATCH'] = '/tmp'
    try:
        os.environ['SLURM_SUBMIT_DIR']
    except KeyError:
        os.environ['SLURM_SUBMIT_DIR'] = '/home/qgm/grendel'

    stepdir = os.path.join(os.environ['LOCALSCRATCH'],
                           'step{}'.format(step))
    lpdir = os.path.join(stepdir, 'prot')
    shutil.copytree(pdir, lpdir)
    pfiles = glob.glob('{}/*.txt'.format(lpdir))
    opt = Option.Option(os.path.join(cfg['optsdir'],
                                     'step{}_opts'.format(step)))

    # Build dict of {(protid, lineno): pattern}
    bonds = {}
    for pfile in pfiles:
        pid = os.path.splitext(os.path.basename(pfile))[0]
        tfile = os.path.join(cfg['tertdir'], '{}.txt'.format(pid))
        if not os.path.exists(tfile) or cfg['max_tbond_level'] == -1:
            tfile = None
        prot = Protein.Protein(pid)
        prot.fromfiles(pfile, tfile)
        for bnd in prot.hbonds:
            bonds[(pid, bnd.lineno)] = str(prot.findpattern2(bnd, opt))

    # Now load assess file
    with open(af) as f:
        ls = f.readlines()
    lines = [l for l in ls
             if (not l.startswith('#')) and (len(l.strip()) > 0)]
    scores = {}
    for line in lines:
        pat, m, e, s, _, _, _, _, _, _ = line.split()
        mode = map(int, m.split(','))
        n = float(e.split(';')[0])
        v = map(float, e.split(';')[1].split(','))
        ev = [n] + v
        score = float(s)
        scores[pat] = (mode, ev, score)

    res = {}
    for bond in bonds:
        if bonds[bond] in scores:
            res[bond] = (bonds[bond],) + scores[bonds[bond]]

    outf = os.path.join(os.environ['SCRATCH'], 'step{}_score.pkl'.format(step))
    with open(outf, 'wb') as o:
        cPickle.dump(res, o, -1)


def listofbonds(pdir):
    res = []
    for f in glob.glob('{}/*.txt'.format(pdir)):
        p, _ = os.path.splitext(os.path.basename(f))
        with open(f) as fh:
            for lineno, line in enumerate(fh):
                cols = line.split()
                if not re.search("__[US][US]$",
                                 cols[cfg['hbond']['flags_col']]):
                    continue
                res.append((p, lineno+1))
    return res


def main(pdir, mstep, outf):
    # Find local pattern around each bond in each test protein for each opts
    # files, and get the results together with the matching line
    # from the corresponding _assess file.
    afs = glob.glob(afpattern)
    for af in afs:
        step = int(os.path.basename(af).split('_')[0].lstrip('step'))
        if os.path.exists(
                os.path.join(
                    os.environ['SCRATCH'], 'step{}_score.pkl'.format(step))):
            print('step{}_score file exists. Skipping.'.format(step))
            continue
        if step > mstep:
            continue
        job_server.submit(runstep,
                          args=(step, pdir, af),
                          group='step{}'.format(step))

    remaining = listofbonds(pdir)

    job_server.wait()

    filestodel = []
    preds = {}
    for step in range(mstep, -1, -1):
        resf = os.path.join(os.environ['SCRATCH'],
                            'step{}_score.pkl'.format(step))
        with open(resf, 'rb') as fh:
            res = cPickle.load(fh)
        with open(resf + '.read', 'w') as w:
            w.write('read.')
        newrem = []
        for bond in remaining:
            try:
                newval = res[bond] + (step,)
            except KeyError:
                # Can't find the bond in clustering result. Search
                # again in the next step (which is step-1).
                newrem.append(bond)
                continue
            if newval[3] == 100:
                # We've found a singleton cluster. No further search
                # needed for this bond.
                preds[bond] = newval
                continue
            try:
                oldval = preds[bond]
            except KeyError:
                # This bond has not been predicted before. Add this
                # to predictions but keep searching.
                preds[bond] = newval
                newrem.append(bond)
                continue
            if newval[3] > oldval[3]:
                # This bond has been predicted before, but the new
                # prediction is better. Add this to predictions but
                # keep searching.
                preds[bond] = newval
            newrem.append(bond)

        # Update list of remaining bonds before the next step
        if len(newrem) == 0:  # We are done with predictions
            break
        remaining = newrem

        filestodel.append(resf)

        with open(resf + '.done', 'w') as w:
            w.write('processed.')

    # Now we have a dict of predictions. Compute the diff between
    # the predicted and actual rotations, and write the results.
    output = []
    for f in glob.glob('{}/*.txt'.format(pdir)):
        with open(f) as fh:
            lines = fh.readlines()
        for lineno, line in enumerate(lines):
            try:
                prot, _ = os.path.splitext(os.path.basename(f))
                pred = preds[(prot, lineno+1)]
            except KeyError:
                # Out list of predictions did not contain the
                # particular bond, probably because local pattern
                # could not be determined. We do the best we can.
                continue
            phi, x, y, z = pred[2]
            # Make sure (x,y,z) is a unit vector.
            n = math.sqrt(x**2 + y**2 + z**2)
            x, y, z = x/n, y/n, z/n
            guessrot = conv.ev2mat(phi, [x, y, z])

            x2, y2, z2, phi2 = map(float, line.split()[15:19])
            # Make sure (x,y,z) is a unit vector.
            n2 = math.sqrt(x2**2 + y2**2 + z2**2)
            x2, y2, z2 = x2/n2, y2/n2, z2/n2
            truerot = conv.ev2mat(phi2, [x2, y2, z2])

            diff = conv.so3dist(truerot, guessrot)
            item = ('{pr}\t{l}\t{s:.4f}\t{d:.6f}\t'
                    '{x:.6f},{y:.6f},{z:.6f}\t'
                    '{pt}\t{st}').format(
                        pr=prot, l=lineno+1, s=pred[3], d=diff,
                        x=x*phi, y=y*phi, z=z*phi,
                        pt=pred[0], st=pred[4]
                    )
            output.append(item)

    with open(outf, 'w') as f:
        f.write('\n'.join(output) + '\n')

    # All done, remove temp files.
    for ftodel in filestodel:
        os.remove(ftodel)
        os.remove(ftodel + '.done')
        os.remove(ftodel + '.read')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protdir', help='Directory containing '
                        'test protein files')
    parser.add_argument('maxstep', type=int,
                        help='Maximum step number '
                        'to consider when making predictions.')
    parser.add_argument('outfile', help='Output file')
    args = parser.parse_args()
    main(args.protdir, args.maxstep, args.outfile)
