#!/usr/bin/env python
#
# File: Option.py
#
# Description: Container class for flp options.
#
# Author: Yuki Koyanagi
#
from collections import namedtuple


class Option(namedtuple('Option', 'window remotes twists residue')):
    """
    Container class for flp options. Load options from _opt file.
    """
    def __new__(cls, optfile):
        with open(optfile) as fh:
            for line in fh:
                try:
                    o, v = line.split()
                except ValueError:
                    continue
                if o == '--window-size':
                    win = int(v)
                elif o == '--nearby-remotes':
                    rem = int(v)
                elif o == '--nearby-twists':
                    twi = int(v)
                elif o == '--residue-scheme':
                    res = int(v)
        self = super(Option, cls).__new__(cls, win, rem, twi, res)
        return self
