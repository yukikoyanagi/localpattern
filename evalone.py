#!/usr/bin/env python
#
# File: eval.py
#
# Description: A fork of eval.py to run single step.
#
# Author: Yuki Koyanagi
# History:
#
import os
import subprocess
import math
import cPickle as pickle
import conv
import argparse


def appendsingles(step, workdir):
    f = os.path.join(workdir, 'singles.pkl')
    with open(f, 'rb') as s:
        singles = pickle.load(s)

    if len(singles) == 0:
        return

    lines = []
    for pattern in singles:
        lines.append('# Processing step{s}/{p}_Summary2.txt'.format(
            s=step, p=pattern
        ))
        # lines.append('# pattern {p}, length {l}'.format(
        #     p='_'.join(pattern.split('_')[0:2]), l=pattern.split('_')[-1]
        # ))
        lines.append('# total observations: 1')
        lines.append('# #clusters: 1')
        lines.append('#\tid\t#obs\tpeakness\tmp-ratio\tmode')
        x, y, z, phi = map(float, singles[pattern])
        # Make sure (x,y,z) is a unit vector.
        # These are from raw data, which may have rounding error
        n = math.sqrt(x**2 + y**2 + z**2)
        x, y, z = x/n, y/n, z/n
        bx, by, bz = conv.ev2box(phi, [x, y, z])
        lines.append('#\t1\t1\t{:f}\t{:f}\t{x},{y},{z}'.format(
            0, 0, x=bx, y=by, z=bz
        ))
        lines.append('#')
        lines.append('# Assessment: Single observation.')
        lines.append('{pat}\t{ix},{iy},{iz}\t'
                     '{p};{fx},{fy},{fz}\t{s:.4f}\t'
                     '{a:.2f}\t{b:.2f}\t{c:.2f}\t'
                     '{d:.2f}\t{e:.2f}\t{f:.2f}'.format(
                         pat=pattern,
                         ix=bx, iy=by, iz=bz,
                         p=phi, fx=x, fy=y, fz=z, s=100,
                         a=0, b=0, c=0, d=0, e=0, f=0
                     ))
        lines.append('#')

    assessf = os.path.join(workdir, 'step{}_assess'.format(step))
    with open(assessf, 'a') as f:
        f.write('\n'.join(lines) + '\n')

def runeval(sumf, step, perlscript, optsdir):
    with open(sumf, 'rb') as f:
        summary = pickle.load(f)

    lines = []
    for pattern in summary:
        prep = 'step{s}/{p}_Summary2.txt\t'.format(
            s=step,
            p=pattern
        )
        # Exclude the last item which is empty when we split at
        # '\n' character
        lst = summary[pattern].split('\n')[:-1]
        for l in lst:
            lines.append(prep + l)

    tempf = 'step{}_summary2.txt'.format(step)
    with open(tempf, 'w') as o:
        o.write('\n'.join(lines))

    try:
        subprocess.check_output([perlscript,
                                 tempf,
                                 '--flp-out-dir', optsdir,
                                 '--output-dir', os.getcwd()],
                                stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        print 'Error while processing {}'.format(tempf)
        print cpe.output
        raise

    outdir = os.path.dirname(sumf)
    if os.path.exists(os.path.join(outdir, 'singles.pkl')):
        appendsingles(step, outdir)
    
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('summaryfile',
                        help='summary.pkl file')
    parser.add_argument('step', type=int,
                        help='Step number')
    parser.add_argument('perlscript',
                        help='eval.pl file')
    parser.add_argument('optsdir',
                        help='Directory with _opts files')
    args = parser.parse_args()

    runeval(args.summaryfile, args.step, args.perlscript, args.optsdir)
