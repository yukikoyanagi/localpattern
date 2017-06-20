#!/usr/bin/env python
#
# File: conv.py
#
# Time-stamp: <2017-03-03 10:36:31 au447708>
#
# Description: Helper functions to perform Euler vector <-> rotation
#  matrix conversions.
#  TODO: Accept non-unit vector as well as angle & unit vector combi.
#
# Author: Yuki Koyanagi
# History:
#
import math
import numpy as np

pi = math.pi


def coord2idx(c):
    '''Box coordinates -> box id's converrion'''
    # The boxes are numbered 1..81, with box i covering [-pi +
    # 2pi*(i-1)/81, -pi + 2pi*i/81]. Hence the index corresponding to
    # a given coordinate w is ceil((w+pi)*81/2pi).
    return map(lambda x: int(math.ceil((x + pi)*81 / (2*pi))), c)


def idx2coord(c):
    '''Box id's -> box coorindates conversion'''
    return map(lambda x: -pi + 2*pi*(x - 0.5)/81, c)


def so3dist(A, B):
    '''Returns angle between two rotation matrices'''
    t = np.dot(A.transpose(), B).trace()
    t = max(-1.0, min(t, 3.0))  # make sure -1<=t<=3
    return math.acos((t - 1)/2)


def ev2mat(angle, u):
    """
    Convert angle, unit vector pair (u is an iterable of x,y,z coords)
    to the corresponding rotation matrix.
    """
    l = math.sqrt(sum(map(lambda x: x**2, u)))
    u = map(lambda x: x/l, u)
    k = np.array([[0, -u[2], u[1]],
                  [u[2], 0, -u[0]],
                  [-u[1], u[0], 0]])
    k2 = np.dot(k, k)
    return np.identity(3) + math.sin(angle)*k + (1-math.cos(angle))*k2


def mat2ev(m):
    '''
    Convert a rotation matrix to the corresponding angle, unit
    vector pair
    '''
    tr = np.trace(m)
    angle = math.acos((tr - 1)/2)
    if angle == 0:
        return 0, (1, 0, 0)

    if math.sin(angle) > 0.0001:
        u = np.array([m[2, 1]-m[1, 2],
                      m[0, 2]-m[2, 0],
                      m[1, 0]-m[0, 1]])
        n = math.sqrt(np.dot(u, u))
        # The original code prints some warning to stderr if the
        # next condition is met, and then ignores it. So don't
        # throw an exception here.
        # if abs(n-2*math.sin(angle)) > 0.001:
        #     raise Exception
        return angle, u/n

    if tr < 0:
        B = (m + np.identity(3)) / 2
        u = np.array(map(math.sqrt, B.diagonal()))
        if u[0] > 0:
            u[1] = u[1] * math.copysign(1, m[0, 1])
            u[2] = u[2] * math.copysign(1, m[0, 2])
        else:
            u[2] = u[2] * math.copysign(1, m[1, 2])

        return angle, u/math.sqrt(np.dot(u, u))


def makemagicmat():
    c = [32, 70, 31]
    u = idx2coord(c)
    n = math.sqrt(np.dot(u, u))
    return np.transpose(ev2mat(n, map(lambda x: x/n, u)))


def box2ev(c):
    u = idx2coord(c)
    angle = math.sqrt(np.dot(u, u))
    M = np.dot(np.transpose(makemagicmat()),
               ev2mat(angle, map(lambda x: x/angle, u)))
    return mat2ev(M)


def ev2box(angle, u):
    B = ev2mat(angle, u)
    M = np.dot(makemagicmat(), B)
    a, v = mat2ev(M)
    return coord2idx(map(lambda x: x*a, v))
