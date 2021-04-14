#!/usr/bin/env python
#
# File: gen_patrot.py
#
# Description: Generate local pattern & rotation for specified protein
# and step numbers. The output is saved in {prot_id}_{step#}.patrot files
# in the specified output directory.
#
# Author: Yuki Koyanagi
#
import os
import argparse
import os.path
from collections import namedtuple

from config import cfg
import Protein
import Option


def main(pf, steps, od):
    pid = os.path.splitext(os.path.basename(pf))[0]
    if '-' in steps:
        s, e = steps.split('-')
        s = int(s)
        e = int(e)
    else:
        s = int(steps)
        e = s
        
    outfile = os.path.join(od, '{}.patrot'.format(pid))
    res = []
    
    for step in range(s, e+1):
        of = os.path.join(cfg['optsdir'], 'step{}_opts'.format(step))
        opt = Option.Option(of)
        tertf = os.path.join(cfg['tertdir'], pid + '.txt')
        if not os.path.exists(tertf):
            # Corresponding tertiary file may not be present...
            tertf = None
        prot = Protein.Protein(pid)
        prot.fromfiles(pf, tertf)

        for bond in prot.hbonds:
            pat = str(prot.findpattern2(bond, opt))
            res.append('\t'.join([str(step), pat, '{},{},{}:{}'.format(*bond.rotation)]))
            
        del prot

    with open(outfile, 'w') as fh:
        fh.write('\n'.join(res))

        
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prot_file',
                        help='Protein hbond file.')
    parser.add_argument('steps',
                        help='Steps (parameter combinations) to '
                        'generate patterns for. A single step '
                        'number or a range of steps (n-m) can '
                        'be specified.')
    parser.add_argument('out_dir',
                        help='Output directory.')
    args = parser.parse_args()
    main(args.prot_file, args.steps, args.out_dir)

