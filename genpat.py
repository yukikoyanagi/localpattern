#!/usr/bin/env python

import os.path
import cPickle
import argparse

import Protein
import Option

def main(protdir, lstf, optfile, outfile):
    pfiles = []
    with open(lstf) as fh:
        for line in fh:
            pfiles.append('{}/{}'.format(os.path.normpath(protdir),
                                         line.strip()))
    opt = Option.Option(optfile)

    # Build dict of {(protid, lineno): pattern}
    bonds = {}
    for pfile in pfiles:
        pid = os.path.splitext(os.path.basename(pfile))[0]
        prot = Protein.Protein(pid)
        prot.fromfiles(pfile, None)
        for bnd in prot.hbonds:
            bonds[(pid, bnd.lineno)] = str(prot.findpattern2(bnd, opt))

    with open(outfile, 'wb') as fh:
        cPickle.dump(bonds, fh, -1)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protdir',
                        help='Directory containing protein files')
    parser.add_argument('lstf',
                        help='File with a list of protein files')
    parser.add_argument('optfile',
                        help='Option file (step##_opts)')
    parser.add_argument('outfile',
                        help='Output pkl file')
    args = parser.parse_args()
    main(args.protdir, args.lstf, args.optfile, args.outfile)

