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
from itertools import groupby
from operator import itemgetter
import conv
import Pattern

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
                rot.append(float(cols[cfg['hbond']['rot_phi']]))
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
                rot = u + [phi]
                vdw = float(cols[cfg['tbond']['vdw_col']])
                self.tbonds.append(Tbond(lidx, ridx, res, rot, vdw))

    def findpattern2(self, center, opts):
        cd = center.donor
        ca = center.acceptor
        tw = self.istwisted(center.rotation)
        cl = (cd + 3)/3 - (ca + 1)/3
        if abs(cl) > 6:
            cl = 'L'
        cbond = Pattern.Bond('H', cd, ca, tw, None)
        dseg = range(cd-opts.window, cd+opts.window+1)
        aseg = range(ca-opts.window, ca+opts.window+1)
        rpat = Pattern.Pattern(segments=[dseg, aseg],
                               bonds=[cbond],
                               residue=center.residue,
                               rotation=center.rotation,
                               length=cl)
        rpat.handleresidue(opts.residue)
        hmax = cfg['max_hbond_level']
        tmax = cfg['max_tbond_level']

        if tmax > 0:
            for b, inend, otherend in self.getogtbonds(rpat):
                self.addt(rpat, b, otherend,
                          opts.window - min(map(abs, [inend-cd, inend-ca]))-1,
                          tmax - 1)

        # Now we need to clean up rpat
        # First we remove duplicates
        atoms = sorted(set([atom for seg in rpat.segments for atom in seg]))
        # then we split the sequence whenever there is a gap.
        # Technically the adjacent atoms may not be connected to each other
        # in a pattern; if atom a and b lie adjacent to each other, but both
        #  are window_size away from the central bond, then they should not
        # be connected. But in practice this is not a problem, since two
        # patterns with the same opts parameters cannot differ by this
        # 'missing' edge.
        rpat.segments = [map(itemgetter(1), g)
                         for k, g
                         in groupby(enumerate(atoms), lambda (i, x): i-x)]

        # Add the bonds internal to the pattern
        tbs2add = [b for b in self.tbonds
                   if (b.left in atoms) and (b.right in atoms)]
        for b in tbs2add:
            rpat.bonds.append(Pattern.Bond('T', b.left, b.right,
                                           self.istwisted(b.rotation), b.vdw))
        hbs2add = [b for b in self.hbonds
                   if (b.donor in atoms) and (b.acceptor in atoms)]
        for b in hbs2add:
            rpat.bonds.append(Pattern.Bond('H', b.donor, b.acceptor,
                                           self.istwisted(b.rotation), None))
        # ... and de-duplicate and sort
        ubonds = set(rpat.bonds)
        rpat.bonds = sorted(ubonds, key=itemgetter(1))

        # Apply --nearby-remotes parameter. We only apply this to hbonds.
        # Remote bonds are the bonds, where only one of the ends is inside
        # the pattern.
        remotes = []
        for hbnd in self.hbonds:
            b = Pattern.Bond('H', hbnd.donor, hbnd.acceptor,
                             self.istwisted(hbnd.rotation), None)
            if b not in rpat.bonds:
                if (b.start in atoms
                    and any([
                        rpat.inwindow(c, b.start, opts.remotes,
                                      hlimit=cfg['max_hbond_level'],
                                      tlimit=cfg['max_tbond_level'])
                        for c in [center.donor, center.acceptor]
                            ])):
                    remotes.append(Pattern.Bond(b.type, b.start,
                                                -99, b.twisted, b.vdw))
                elif (b.end in atoms
                      and any([
                        rpat.inwindow(c, b.end, opts.remotes,
                                      hlimit=cfg['max_hbond_level'],
                                      tlimit=cfg['max_tbond_level'])
                        for c in [center.donor, center.acceptor]
                      ])):
                    remotes.append(Pattern.Bond(b.type, -99,
                                                b.end, b.twisted, b.vdw))
        rpat.bonds += remotes

        # Apply --nearby-twists parameter. Note nearby-twists takes only
        # three values; -1, 0 and window-size. So if twists > 0, we can
        # include twist information on all bonds.
        if opts.twists <= 0:
            newbonds = []
            for bond in rpat.bonds:
                if (bond.start == center.donor and
                        bond.end == center.acceptor and
                        opts.twists == 0):
                    newbonds.append(bond)
                else:
                    newbonds.append(Pattern.Bond(bond.type, bond.start,
                                                 bond.end, False, bond.vdw))
            rpat.bonds = newbonds

        # We need to replace atom index by the local index, i.e.
        # 0, 1, 2,... for the donor-segment for the central bond
        # 100, 101, ... for the acceptor-segment for the central bond
        # 200, 201, ... for the segment connected to the '0' segment at the
        # smallest index
        # ...and so on
        cbond = Pattern.Bond('H', center.donor, center.acceptor,
                             self.istwisted(center.rotation), None)
        lpat = rpat.localise(cbond)

        return lpat



    def findpattern(self, center, opts):
        """
        Return a Pattern object around *center* Hbond using *opts*.
        :type center: Hbond
        :param center: hbond object for the central bond of the pattern
        :param opts: Option object
        :return: a Pattern object
        """
        dseg = Pattern.Segment()
        dseg.append(center.donor)
        aseg = Pattern.Segment()
        aseg.append(center.acceptor)

        twist = self.istwisted(center.rotation)
        bnd = Pattern.Bond('H', center.donor, center.acceptor, twist, None)
        rpat = Pattern.Pattern(segments=[dseg, aseg],
                               bonds=[bnd],
                               residue=center.residue,
                               rotation=center.rotation)
        rpat.handleresidue(opts.residue)

        # We grow rpat by opts.window. This ignores config parameters
        # 'max_hbond_level' and 'max_tbond_level'.
        i = opts.window
        while i:
            rpat = self.grow(rpat)
            i -= 1

        # Clean up rpat
        # We need to limit the pattern according to the configuration and
        # _opts parameters
        # First trim the pattern according to the config parameters
        # 'max_hbond_level' and 'max_tbond_level'.
        rpat.trim(Pattern.Bond('H', center.donor, center.acceptor,
                               self.istwisted(center.rotation), None),
                  cfg['max_hbond_level'],
                  cfg['max_tbond_level'],
                  opts)

        # Apply --nearby-remotes parameter. We only apply this to hbonds.
        # Remote bonds are the bonds, where only one of the ends is inside
        # the pattern.
        atoms = [a for seg in rpat.segments for a in seg]
        remotes = []
        for hbnd in self.hbonds:
            b = Pattern.Bond('H', hbnd.donor, hbnd.acceptor,
                             self.istwisted(hbnd.rotation), None)
            if b not in rpat.bonds:
                if (b.start in atoms
                    and any([
                        rpat.inwindow(c, b.start, opts.remotes,
                                      hlimit=cfg['max_hbond_level'],
                                      tlimit=cfg['max_tbond_level'])
                        for c in [center.donor, center.acceptor]
                            ])):
                    remotes.append(Pattern.Bond(b.type, b.start,
                                                -99, b.twisted, b.vdw))
                elif (b.end in atoms
                      and any([
                        rpat.inwindow(c, b.end, opts.remotes,
                                      hlimit=cfg['max_hbond_level'],
                                      tlimit=cfg['max_tbond_level'])
                        for c in [center.donor, center.acceptor]
                      ])):
                    remotes.append(Pattern.Bond(b.type, -99,
                                                b.end, b.twisted, b.vdw))
        rpat.bonds += remotes

        # Apply --nearby-twists parameter. Note nearby-twists takes only
        # three values; -1, 0 and window-size. So if twists > 0, we can
        # include twist information on all bonds.
        if opts.twists <= 0:
            newbonds = []
            for bond in rpat.bonds:
                if (bond.start == center.donor and
                        bond.end == center.acceptor and
                        opts.twists == 0):
                    newbonds.append(bond)
                else:
                    newbonds.append(Pattern.Bond(bond.type, bond.start,
                                                 bond.end, False, bond.vdw))
            rpat.bonds = newbonds

        # We need to replace atom index by the local index, i.e.
        # 0, 1, 2,... for the donor-segment for the central bond
        # 100, 101, ... for the acceptor-segment for the central bond
        # 200, 201, ... for the segment connected to the '0' segment at the
        # smallest index
        # ...and so on
        cbond = Pattern.Bond('H', center.donor, center.acceptor,
                             self.istwisted(center.rotation), None)
        lpat = rpat.localise(cbond)

        return lpat

    def istwisted(self, rotation):
        """
        We can determine whether a bond is twisted by looking at the last
        entry (i.e. (3,3) entry).
        :param rotation: a list of 9 entries in SO3 rotation matrix
        :return: Boolean
        """
        if len(rotation) == 9:
            # We assume it's an SO3 matrix
            return rotation[-1] < 0
        elif len(rotation) == 4:
            # We assume it's of the form (x, y, z, rot)
            l = conv.mat2lst(conv.ev2mat(rotation[3], tuple(rotation[:3])))
            return l[-1] < 0

    def getbondsat(self, pos):
        """
        Get a list of bonds (Hbond or Tbond) at the *pos* position.
        :param pos: int; atom id
        :return: a list of bonds
        """
        result = []
        for hbnd in self.hbonds:
            if hbnd.donor == pos or hbnd.acceptor == pos:
                result.append(hbnd)
        for tbnd in self.tbonds:
            if tbnd.left == pos or tbnd.right == pos:
                result.append(tbnd)
        return result

    def grow(self, pat):
        """
        Grow a Pattern object by one; i.e. return a Pattern object which
        include the edge/atom pairs which are adjacent to the given Pattern
        :param pat: Pattern object
        :return: A 'grown' Pattern object
        """
        newsegs = pat.segments[:]
        newbonds = pat.bonds[:]
        for seg in pat.segments:
            for edge in [seg[0], seg[-1]]:
                bonds = self.getbondsat(edge)
                for bond in bonds:
                    if isinstance(bond, Hbond):
                        b = Pattern.Bond('H', bond.donor, bond.acceptor,
                                         self.istwisted(bond.rotation), None)
                    elif isinstance(bond, Tbond):
                        b = Pattern.Bond('T', bond.left, bond.right,
                                         self.istwisted(bond.rotation),
                                         bond.vdw)
                    if b not in newbonds:
                        newbonds.append(b)
                        # We simply add both ends as new segments to the
                        # pattern, and clean up later.
                        newsegs.append([b[1]])
                        newsegs.append([b[2]])
                if seg.index(edge) == 0:  # this is a left edge
                    newsegs.append([edge - 1] + seg)
                if seg.index(edge) == len(seg) - 1:  # this is a right edge
                    newsegs.append(seg + [edge + 1])

        # Now we need to clean up the pattern.
        # let's first remove duplicate atoms
        atoms = set([atom for segment in newsegs for atom in segment])
        seq = sorted(atoms)
        # then we split the sequence whenever there is a gap.
        # Technically the adjacent atoms may not be connected to each other
        # in a pattern; if atom a and b lie adjacent to each other, but both
        #  are window_size away from the central bond, then they should not
        # be connected. But in practice this is not a problem, since two
        # patterns with the same opts parameters cannot differ by this
        # 'missing' edge.
        pat.segments = [map(itemgetter(1), g)
                        for k, g
                        in groupby(enumerate(seq), lambda (i, x): i-x)]

        # We need to add the bonds which have both start and end atoms at
        # the ends of segments. They do not get captured by the above,
        # but are clearly internal to the pattern.
        ends = [p for seg in pat.segments for p in [seg[0], seg[-1]]]
        for end in ends:
            for bond in self.getbondsat(end):
                if isinstance(bond, Hbond):
                    b = Pattern.Bond('H', bond.donor, bond.acceptor,
                                     self.istwisted(bond.rotation), None)
                elif isinstance(bond, Tbond):
                    b = Pattern.Bond('T', bond.left, bond.right,
                                     self.istwisted(bond.rotation),
                                     bond.vdw)

                if all([b[1] in ends,
                        b[2] in ends,
                        b not in newbonds]):
                    newbonds.append(b)
        pat.bonds = sorted(newbonds, key=itemgetter(1))

        return pat

    def getogtbonds(self, pat):
        """
        Return a list of 'outgoing' tbonds in *pat*. 'Outgoing' tbonds are
        the tbonds which have only one of the ends in the pattern.
        :return: List of tuple containing;
         bond: Pattern.Bond instance for an outgoing tbond
         inend: int; index of the end which is in the pattern
         otherend: int; index of the end outside the pattern
        """
        atoms = [a for seg in pat.segments for a in seg]
        results = []
        for tb in self.tbonds:
            if (tb.left in atoms) is not (tb.right in atoms):  # XOR
                bnd = Pattern.Bond('T', tb.left, tb.right,
                                   self.istwisted(tb.rotation), tb.vdw)
                # i is the end inside the pattern, since True evaluates to 1
                # and False to 0
                i = tb.left * (tb.left in atoms) +\
                    tb.right * (tb.right in atoms)
                o = tb.left + tb.right - i
                results.append((bnd, i, o))
        return results

    def addt(self, pat, bond, c, win, tlevel):
        """
        'Grow' *pat* by adding T-section(s) to it. T-section is a pairing of
        *bond* and a segment on the other side of the *bond*. *win* and
        *tlevel* provide constraints.
        :param pat: Pattern.Pattern object to be grown
        :param bond: Pattern.Bond instance; outgoing t-bond
        :param c: ind; the other end of *bond*
        :param win: int; Remaining window size after adding the bond
        :param tlevel: int; Remaining max_tbond_level after adding the bond
        :return: expanded Pattern.Pattern instance
        """
        if win >= 0:
            pat.bonds.append(bond)
            newseg = Pattern.Segment(range(c-win, c+win+1))
            pat.segments.append(newseg)
        if tlevel > 0 and win > 0:
            for b, i, o in self.getogtbonds(pat):
                if i in newseg:
                    self.addt(pat, b, o, win-abs(i-c)-1, tlevel-1)
        return pat
