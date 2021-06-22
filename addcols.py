#!/usr/bin/env python
#
# File: addcols.py
#
# Time-stamp: <2017-02-23 11:01:23 au447708>
#
# Description: Adds required cols to raw protein file.
#
# Author: Yuki Koyanagi
# History:
#
import os
import sys
import argparse


def main(protf, codepath=None, outdir=None):
    if codepath:
        cpath = os.path.abspath(codepath)
        sys.path.append(cpath)
    import conv

    with open(protf) as fh:
        data = fh.readlines()
    bonds = []
    for line in data:
        bond = line.split()
        bond = bond + ['?']*5
        x, y, z, phi = map(float, bond[15:19])
        M = conv.ev2mat(phi, [x, y, z])
        bond = bond + map('{:f}'.format, sum(M.tolist(), []))
        bonds.append('\t'.join(bond))

    if outdir:
        protf = os.path.join(outdir, os.path.basename(protf))
    with open(protf, 'w') as fh:
        fh.write('\n'.join(bonds))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protf')
    parser.add_argument('-p', '--codepath', type=str, default=None)
    parser.add_argument('-o', '--outdir', default=None)
    args = parser.parse_args()
    main(args.protf, args.codepath, args.outdir)
