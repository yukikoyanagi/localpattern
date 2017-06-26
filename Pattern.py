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

"""Bond represents a h- or t-bond.
type: either h or t
start, end: id of start (donor/left) and end (acceptor/right) atoms
ex. 6: 7th atom in the first segment, 102: 3rd atom in the second segment
twisted: Boolean, whether the bond is twisted
"""
Bond = namedtuple('Bond', 'type start end twisted')

"""Atom represents an atom along a backbone segment.
position: the atom's position along the backbone segment
residue: the residue class (where appropriate) for the atom
"""
Atom = namedtuple('Atom', 'position residue')


class Segment(list):
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
        rotation: the rotation associated to the central bond
    """
    def __init__(self, segments=None, bonds=None, rotation=None):
        if not segments:
            segments = []
        self.segments = segments
        if not bonds:
            bonds = []
        self.bonds = bonds
        self.rotation = rotation

    def __eq__(self, other):
        return self.segments == other.segments and self.bonds == other.bonds

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.segments, self.bonds))
