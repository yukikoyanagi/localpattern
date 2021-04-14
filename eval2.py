#!/usr/bin/env python
#
# File: eval2.py
#
# Description:
#

import os, subprocess, math, cPickle, argparse, os.path

from config import cfg

def main(inf, dbg):
    minn = cfg['min_for_cluster']
    outd = os.path.dirname(inf)
    step = int(os.path.basename(inf).split('_')[0][4:])
    
    with open(inf, 'rb') as fh:
        summary = cPickle.load(fh)

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

    tempf = os.path.join(outd, 'step{}_summary2.txt'.format(step))
    with open(tempf, 'w') as fh:
        fh.write('\n'.join(lines))

    optsd = cfg['optsdir']
    perlscript = 'eval.pl'

    try:
        subprocess.check_output(['perl', perlscript,
                                 tempf,
                                 '--flp-out-dir', optsd,
                                 '--output-dir', outd],
                                stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        print 'Error while processing {}'.format(tempf)
        print cpe.output

    if not dbg:
        os.remove(tempf)


    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file',
                        help='Input pickle file, step##_summary.pkl.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Keep all intermediate files.')
    args = parser.parse_args()
    main(args.in_file, args.debug)
