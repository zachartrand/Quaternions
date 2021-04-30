# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 21:23:28 2021

@author: Zach Chartrand <zachartrand999@gmail.com>

Similar to the built-in module cmath, this module has definitions of
mathematical functions expanded to work with quaternions.
"""

__all__ = ['exp', 'log', 'pow', 'rotate3d']

from math import exp as _exp, cos as _cos, sin as _sin, log as _log
from math import log1p as _log1p, acos as _acos, pow as _pow

from typing import Iterable

from quaternions import Quaternion

def exp(q: Quaternion) -> Quaternion:
    """Returns the exponential of a quaternion."""
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_exp(a), 0, 0, 0)
        else:
            v = q.vector
            return _exp(a)*(_cos(abs(v)) + v/abs(v) * _sin(abs(v)))

def log(q: Quaternion) -> Quaternion:
    """Return the natural logarithm of a quaternion."""
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            if abs(a - 1) < 0.1:
                return Quaternion(_log1p(a - 1), 0, 0, 0)
            else:
                return Quaternion(_log(a), 0, 0, 0)
        else:
            v = q.vector
            if abs(abs(q) - 1) < 0.1:
                return _log1p(abs(q) - 1) + v/abs(v)*_acos(a/abs(q))
            else:
                return _log(abs(q)) + v/abs(v)*_acos(a/abs(q))

def pow(q: Quaternion, x: int or float) -> Quaternion:
    """Return q**x, where x is a real exponent."""
    if (isinstance(q, Quaternion) and (isinstance(x, int)
            or isinstance(x, float))):
        if x // 1 == x:
            for _ in range(x - 1):
                q *= q
            return q
        else:
            a = q.scalar
            if q.is_scalar():
                return Quaternion(_pow(a, x), 0, 0, 0)
            else:
                v = q.vector
                theta = _acos(a/abs(q))
                return _pow(abs(q), x) * exp(v.unit_vector()*x*theta)

def rotate3d(point: tuple or list, angle: int or float,
             axis: tuple or int=(0, 0, 1), rounding: int=14) -> tuple:
    """
    Takes a point in 3d space represented as a tuple or list of three
    (3) values and rotates it by an angle around a given axis vector.
    The axis of rotation is the z-axis [0, 0, 1] by default. For the
    point and axis parameters, if only one value is given, the value
    will be assumed to be an x-coordinate with the y- and z-coordinates
    equal to zero (0). If two values are given, they will be assumed to
    be x- and y-coordinates with the z-coordinate equal to zero (0).
    By default, this function rounds all floats in the result to 14
    decimal places. If you wish to round to a different number of
    places, change the rounding argument to the number of decimal
    places you wish to round to. If you don't want any rounding, set
    the parameter to -1.
    """
    if len(point) <= 3 and len(axis) <= 3:
        i, j, k = _makeListLen3(point)
        p = Quaternion(0, i, j, k)

        i2, j2, k2 = _makeListLen3(axis)
        u = Quaternion(0, i2, j2, k2)
        if abs(u) == 0:
            raise ValueError(
                'The axis to rotate around must be a nonzero vector.')
        if abs(u) != 1:
            u = u/abs(u)
        q = _cos(angle/2) + u*_sin(angle/2)
        p_prime = q * p * q.inverse()
        i_prime, j_prime, k_prime = p_prime.get_vector_components()
        if rounding != -1:
            i_prime, j_prime, k_prime = (
                round(i_prime, rounding), round(j_prime, rounding),
                round(k_prime, rounding))

        return [i_prime, j_prime, k_prime]

def _makeListLen3(i: Iterable[int or float]) -> list:
    if not isinstance(i, list):
        i = list(i)
    while len(i) < 3:
        i.append(0.0)
    return i
