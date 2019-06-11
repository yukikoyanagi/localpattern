#!/usr/bin/env python
#
# File: filter.py
#
# Created by: au447708 on 8/15/17
#
# Description: Filter for a single step

import os
import cPickle
import argparse

def main(inf, n, outf):
    with open(inf, 'rb') as fh:
        pr = cPickle.load(fh)
    outrot = {p: pr[p] for p in pr if len(pr[p]) >= n}
    with open(outf, 'wb') as o:
        cPickle.dump(outrot, o, -1)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Input rotation pkl file')
    parser.add_argument('min', type=int,
                        help='Min. number of occurrences')
    parser.add_argument('output', help='Output file')
    args = parser.parse_args()
    main(args.input, args.min, args.output)
