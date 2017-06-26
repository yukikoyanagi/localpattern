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

class Bond(namedtuple('Bond', 'type start end twisted')):
    """
    Class representing an h- or t-bond.
    Attributes:
        type: either h or t
        start, end: id of start (donor/left) and end (acceptor/right) atoms
        twisted: Boolean, whether the bond is twisted
    """
    __slots__ = ()

    def __eq__(self, other):
        return (self.type == other.type and self.twisted == other.twisted)

    def __ne__(self, other):
        return not self == other


class Atom(namedtuple('Atom', 'position residue')):
    """
    Class representing an atom along a backbone segment.
    Attributes:
        position: the atom's position along the (entire) backbone
        residue: the residue class (where appropriate) for the atom
    """
    __slots__ = ()

    def __eq__(self, other):
        return self.residue == other.residue

    def __ne__(self, other):
        return not self == other


class Segment(list):
    pass


class Pattern(object):
    """
    Local pattern description of a hbond.
    Attributes:
        segments: set of Segment objects.
        bonds: list of Bond objects
    """
    def __init__(self, segments=None, bonds=None):
        if not segments:
            segments = set()
        self.segments = segments
        if not bonds:
            bonds = []
        self.bonds = bonds

    def __eq__(self, other):
        if self.segments != other.segments or self.bonds != other.bonds:
            return False
