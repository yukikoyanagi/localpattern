#!/usr/bin/env python
#
# File: predict.py
#
# Description: Run H-bond rotation prediction for a given protein.
#

import argparse, os, glob, math, cPickle, os.path

import conv
import Protein
import Option
from config import cfg

def runstep(step, prot, af):
    '''
    Runs predictions for the given step. Returns a dict of
    {(protein, lineno): value}, value is a tuple containing
    pattern: string
    mode: 3-tuple of mode box id's
    ev: 4-tuple of rotation vector; angle, x, y, z
    score: score
    '''
    opt = Option.Option(os.path.join(cfg['optsdir'],
                                     'step{}_opts'.format(step)))
    pid = prot.name

    # Build dict of {(protid, lineno): pattern}
    bonds = {}
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

    return res


def main(protf, assessd, outf):
    assessd = assessd.rstrip('/')
    outd = os.path.dirname(outf)
    afs = glob.glob('{}/step*_assess'.format(assessd))
    steps = sorted([int(os.path.basename(af).split('_')[0][4:]) for af in afs])
    pid = os.path.splitext(os.path.basename(protf))[0]
    tfile = os.path.join(cfg['tertdir'], '{}.txt'.format(pid))
    if not os.path.exists(tfile) or cfg['max_tbond_level'] == -1:
        tfile = None
    prot = Protein.Protein(pid)
    prot.fromfiles(protf, tfile)

    # Make score file for each step
    for step in steps:
        af = '{}/step{}_assess'.format(assessd, step)
        step_res = runstep(step, prot, af)
        stepout = os.path.join(outd, '{}_{}_score.pkl'.format(pid, step))
        with open(stepout, 'wb') as fh:
            cPickle.dump(step_res, fh, -1)

    # Make dict of predictions
    remaining = [(pid, bnd.lineno) for bnd in prot.hbonds]
    preds = {}
    for step in sorted(steps, reverse=True):
        resf = os.path.join(outd, '{}_{}_score.pkl'.format(pid, step))
        with open(resf, 'rb') as fh:
            res = cPickle.load(fh)
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

        os.remove(resf)

    # Now we have a dict of predictions. Compute the diff between
    # the predicted and actual rotations, and write the results.
    output = []
    for bnd in prot.hbonds:
        try:
            pred = preds[(pid, bnd.lineno)]
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
        
        x2, y2, z2, phi2 = bnd.rotation
        n2 = math.sqrt(x2**2 + y2**2 + z2**2)
        x2, y2, z2 = x2/n2, y2/n2, z2/n2
        truerot = conv.ev2mat(phi2, [x2, y2, z2])
        
        diff = conv.so3dist(truerot, guessrot)
        item = ('{pr}\t{l}\t{s:.4f}\t{d:.6f}\t'
                '{x:.6f},{y:.6f},{z:.6f}\t'
                '{pt}\t{st}').format(
                    pr=prot.name, l=bnd.lineno, s=pred[3], d=diff,
                    x=x*phi, y=y*phi, z=z*phi,
                    pt=pred[0], st=pred[4]
                )
        output.append(item)

    with open(outf, 'w') as f:
        f.write('\n'.join(output) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protfile',
                        help='Protein Hbond file to predict rotations.')
    parser.add_argument('assessdir',
                        help='Directory with _assess files.')
    parser.add_argument('outfile',
                        help='Output file')
    args = parser.parse_args()
    main(args.protfile, args.assessdir, args.outfile)
