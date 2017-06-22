#!/usr/bin/env python
#
# File: Protein.py
#
# Description: Protein container class, containing a list of Hbonds and Tbonds.
#  Find_Pattern to produce Pattern object.
#
# Author: Yuki Koyanagi
#
import re
from config import cfg
from collections import namedtuple
import conv

Hbond = namedtuple('Hbond', 'donor acceptor residue rotation')
Tbond = namedtuple('Tbond', 'left right residue rotation vdw')


class Protein(object):
    """
    A collection of Hbonds and Tbonds in a protein.
    Attributes:
        name: protein name
        hbonds: list of Hbond objects
        tbonds: list of Tbond objects
    Methods:
        fromfiles:
        findpattern: Return a Pattern object describing local pattern.
    """

    def __init__(self, name):
        """
        Return a Protein object with given *name* and empty hbonds and tbonds.
        """
        self.name = name
        self.hbonds = []
        self.tbonds = []

    def fromfiles(self, hfile, tfile=None):
        """
        Populate self.hbonds and self.tbonds using information in *hfile* and
        *tfile*.
        :param hfile: path to the file containing hbonds
        :param tfile: path to the file containing tbonds
        :return:
        """
        with open(hfile) as hf:
            for line in hf:
                cols = line.split()
                # we only consider intra-chain bonds (indicated by __),
                # and the last two characters must be either U (unique)
                # or S (strong)
                if not re.search("__[US][US]$",
                                 cols[cfg['hbond']['flags_col']]):
                    continue

                donidx = (int(cols[cfg['hbond']['donor_col']]) - 1) * 3
                accidx = (int(cols[cfg['hbond']['accpt_col']]) - 1) * 3 + 2
                residue = [cols[i] for i in cfg['hbond']['res_cols']]
                rot = map(float, [cols[i] for i in cfg['hbond']['rot_cols']])
                self.hbonds.append(Hbond(donidx, accidx, residue, rot))

        if not tfile:
            return
        with open(tfile) as tf:
            for line in tf:
                cols = line.split()
                if cols[cfg['tbond']['flags_col']][2:4] != '__':
                    continue
                lidx = (int(cols[cfg['tbond']['left_col']]) - 1) * 3 + 1
                ridx = (int(cols[cfg['tbond']['right_col']]) - 1) * 3 + 1
                # Only consider the case where lidx < ridx; otherwise we'll
                # have two entries for each tbond.
                if not lidx < ridx:
                    continue
                res = cols[cfg['tbond']['res_col']]
                phi = float(cols[cfg['tbond']['rot_phi']])
                u = map(float, [cols[i] for i in cfg['tbond']['rot_cols']])
                rot = conv.mat2lst(conv.ev2mat(phi, u))
                vdw = float(cols[cfg['tbond']['vdw_col']])
                self.tbonds.append(Tbond(lidx, ridx, res, rot, vdw))



    def findpattern(self, center, opts):
        """
        Return a Pattern object around *center* atom using *opts*.
        :param center: index of the central atom for the pattern
        :param opts: Opts object (?)
        :return: a Pattern object
        """