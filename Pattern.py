#!/usr/bin/env python
#
# File: Pattern.py
#
# Description: Class describing local pattern. Produced by running
#  locpat.Protein.findpattern()
#
# Author: Yuki Koyanagi
#
from collections import namedtuple
from operator import itemgetter


class Bond(namedtuple('Bond', 'type start end twisted vdw')):
    """Bond represents a h- or t-bond.
    type: either H or T
    start, end: id of start (donor/left) and end (acceptor/right) atoms
    ex. 6: 7th atom in the first segment, 102: 3rd atom in the second segment
    twisted: Boolean, whether the bond is twisted. This is always False if
    --neaby-twist is -1.
    vdw: VDW distance. None for h-bond.
    'Remote' bond is represented by '-99' at the 'remote' (i.e. corresponding
    segment not included in the pattern) end.
    """
    __slots__ = ()

    def __eq__(self, other):
        return (self.type == other.type and
                self.start == other.start and
                self.end == other.end and
                self.twisted == other.twisted)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.type, self.start, self.end, self.twisted))


"""Atom represents an atom along a backbone segment.
position: the atom's position along the backbone segment
residue: the residue class (where appropriate) for the atom
"""
# We won't be using the residue information for each atom, only the four
# residues surrounding the central bond. So the atom can simply be
# respresented by an integer.
# Atom = namedtuple('Atom', 'position residue')


class Segment(list):
    """
    This is a list of integers (atom id's) representing a backbone segment.
    """
    pass


class Pattern(object):
    """
    Local pattern description of a hbond. This is simply a container for
    segments and bonds.
    Attributes:
        segments: list of Segment objects. Segments are ordered as follows;
        1. donor-side segment of the central hbond
        2. acceptor-side segment of the central hbond
        Then we move along the 1st segment and list connected segments in the
        order of appearance. Then 2nd segment, and so on.
        bonds: list of Bond objects
        residue: 4 residues around the central bond
        rotation: the rotation associated to the central bond
    """
    def __init__(self, segments=None, bonds=None, residue=None, rotation=None):
        if not segments:
            segments = []
        self.segments = segments
        if not bonds:
            bonds = []
        self.bonds = bonds
        self.residue = residue
        self.rotation = rotation

    def __eq__(self, other):
        return (self.segments == other.segments and
                set(self.bonds) == set(other.bonds) and
                self.residue == other.residue)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.segments, self.bonds, self.residue))

    def trim(self, center, hlimit, tlimit, opts):
        """
        Trim pattern around the *center* bond so that the pattern lies within
        *hlimit* number of hbonds and *tlimit* number of tbonds from the
        central bond.
        :param center: self.Bond instance
        :param hlimit: int;
        :param tlimit: int;
        :return: Trimmed pattern object
        """
        newsegs = []
        for seg in self.segments:
            if self.lieswithin(seg, center, hlimit, tlimit):
                newsegs.append(seg)
        self.segments = newsegs

        # The first check using lieswithin does not catch the case,
        # where the other end of bond was close enough, so that the two
        # segments (one of the central segments and the other) have grown
        # together. But no segment should be longer than 2*w+1, where w is
        # the window size.
        newsegs = []
        for seg in self.segments:
            if center.start in seg or center.end in seg:
                nseg = []
                for a in seg:
                    if min([abs(a - center.start),
                            abs(a - center.end)]) <= opts.window:
                        nseg.append(a)
                newsegs.append(sorted(nseg))
            else:
                newsegs.append(seg)
        self.segments = newsegs

        newbonds = []
        atoms = [a for seg in self.segments for a in seg]
        for bond in self.bonds:
            if bond.start in atoms and bond.end in atoms:
                newbonds.append(bond)
        self.bonds = newbonds

        return

    def lieswithin(self, seg, center, hlimit, tlimit):
        """
        Check if the segment *seg* lies within *hlimit* & *tlimit* of
        H-/T-bonds from the central bond *center*
        :param seg:
        :param center:
        :param hlimit:
        :param tlimit:
        :return: Boolean
        """
        if hlimit < 0 or tlimit < 0:
            return False
        elif center.start in seg or center.end in seg:
            return True
        elif hlimit == 0 and tlimit == 0:
            return False
        else:
            # get all segments connected to seg, and check if any of them is
            #  within the range of central bond.
            bonds = [bond for bond in self.bonds
                     if {bond.start, bond.end} & set(seg)]
            segs = []
            for bond in bonds:
                if bond.start in seg:
                    v = bond.end
                else:
                    v = bond.start
                for segment in self.segments:
                    if v in segment:
                        if bond.type == 'H':
                            hlimit = hlimit - 1
                        else:
                            tlimit = tlimit - 1

                        segs.append((segment, hlimit, tlimit))

            return any([self.lieswithin(s[0], center, s[1], s[2])
                        for s in segs])

    def dist(self, start, end, visited=None):
        """
        Compute the distance between two atoms in a pattern.
        :param start: index of 'from' atom
        :param end: index of 'to' atom
        :param visited: a list of bonds already traversed
        :return: int;
        """
        if visited is None:
            visited = []

        # Select the unique segment that contains start
        seg, = (s for s in self.segments if start in s)
        if end in seg:
            return abs(start - end)
        else:
            bonds = [bond for bond in self.bonds
                     if {bond.start, bond.end} & set(seg) and
                     bond not in visited]
            if len(bonds) == 0:  # We are at a dead end...
                return 999
            l = []
            for bond in bonds:
                if bond.start in seg:
                    s = bond.start
                    e = bond.end
                else:
                    s = bond.end
                    e = bond.start
                l.append((bond, s, e))
            # for each of the connecting bonds, we return
            # (dist from start to the end of bond + 1 (for the bond)
            # + dist from the other end of bond to end) and return the min
            # If all dist is > 999 (i.e. end is unreachable from start),
            # return 999.
            return min([999,
                        min([abs(start - b[1]) +
                             1 +
                             self.dist(b[2], end, visited + [b[0]])
                             for b in l])
                        ])

    def localise(self, center):
        """
        Localise the pattern around the *center* bond. After being
        localised, the bonds will have atom indices given by;
        ex. 6: 7th atom in the 1st segment, 102: 3rd atom in the 2nd segment
        :param center: Bond; central hbond
        :return: localised Pattern object
        """
        # Order the segments & bonds
        ordsegs = {}
        sseg, = (s for s in self.segments if center.start in s)
        eseg, = (s for s in self.segments if center.end in s)
        ordsegs[0] = sseg
        if eseg != sseg:
            ordsegs[1] = eseg

        i = len(ordsegs)
        while len(ordsegs) < len(self.segments):
            newsegs = {}
            for k in ordsegs:
                for j in ordsegs[k]:
                    bs = sorted(
                        [b for b in self.bonds if b.start == j or b.end == j],
                        key=itemgetter(4))
                    for b in bs:
                        if b.start == j:
                            v = b.end
                        else:
                            v = b.start
                        if v != -99:
                            seg, = (s for s in self.segments if v in s)
                            if seg not in ordsegs.values() + newsegs.values():
                                newsegs[i] = seg
                                i = i + 1
            ordsegs.update(newsegs)

        lpat = Pattern()

        for k in ordsegs:
            lpat.segments.append(
                Segment(range(k*100, k*100 + len(ordsegs[k]))))
        for bnd in self.bonds:
            if bnd.start != -99:
                i, = (k for k in ordsegs if bnd.start in ordsegs[k])
                locstart = i*100 + abs(bnd.start - ordsegs[i][0])
            else:
                locstart = bnd.start  # == -99...
            if bnd.end != -99:
                j, = (k for k in ordsegs if bnd.end in ordsegs[k])
                locend = j*100 + abs(bnd.end - ordsegs[j][0])
            else:
                locend = bnd.end
            lpat.bonds.append(Bond(
                bnd.type, locstart, locend, bnd.twisted, bnd.vdw
            ))

        lpat.residue = self.residue
        lpat.rotation = self.rotation
        return lpat

    def handleresidue(self, scheme):
        """
        A residue scheme is a mapping defining groups of "similar" residues.
        The trivial scheme is the one where each residue is mapped to its
        own group. We represent a group simply by its first member.
        It is up to the user to keep track of which scheme was used.
        :param scheme: int;
        :return: classified residues
        """
        schemes = ['XLVIFMAGSCEKRDTYNQHWP',
                   list('XLVIFMAGSCEKRDTYNQHWP'),
                   ["LVIFM", "AGSC", "EKRDTYNQHW", "P", "X"],
                   ["LVIFM", "AGSC", "EKRDTYNQHWP", "X"],
                   ["LVIFMAGSC", "EKRDTYNQHWP", "X"]]
        done = ''
        for res in self.residue:
            if res not in schemes[0]:
                res = 'X'
            newres, = (s[0] for s in schemes[scheme] if res in s)
            done += newres

        self.residue = done
        return
