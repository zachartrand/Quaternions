# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 21:23:28 2021

@author: Zach Chartrand <zachartrand999@gmail.com>

Similar to the built-in module cmath, this module has definitions of
mathematical functions expanded to work with quaternions.
"""

__all__ = ['exp', 'log', 'pow', 'rotate3d']

from math import exp as _exp, cos as _cos, sin as _sin, log as _log
from math import log1p as _log1p, pow as _pow, expm1 as _expm1

from typing import Iterable

from quaternions import Quaternion

def exp(q: Quaternion or int or float) -> Quaternion:
    """Returns the exponential of a quaternion."""
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_exp(a), 0, 0, 0)

        v = q.vector
        return _exp(a)*(_cos(abs(v)) + q.unit_vector()*_sin(abs(v)))

    elif isinstance(q, (int, float)):
        return Quaternion(_exp(q), 0, 0, 0)

def expm1(q: Quaternion) -> Quaternion:
    """
    Return exp(q)-1.

    This function avoids the loss of precision involved in the direct
    evaluation of exp(q)-1 for q with a small norm.
    """
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_expm1(a), 0, 0, 0)

        v = q.vector
        return _expm1(a)*(cos(abs(v)) + q.unit_vector()*_sin(abs(v)))

    elif isinstance(q, (int, float)):
        return Quaternion(_expm1(q), 0, 0, 0)

    return NotImplemented

def log(q: Quaternion or int or float) -> Quaternion:
    """Return the natural logarithm of a quaternion."""
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_log(a), 0, 0, 0)

        angle = q.angle
        return _log(abs(q)) + q.unit_vector()*angle

    elif isinstance(q, (int, float)):
        return Quaternion(_log(q), 0, 0, 0)

    return NotImplemented

def log1p(q: Quaternion or int or float) -> Quaternion:
    """
    Return the natural logarithm of 1+q.

    The result is computed in a way which is accurate for q with a norm
    near zero.
    """
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_log1p(a), 0, 0, 0)

        angle = q.angle
        return _log1p(abs(q)) + q.unit_vector()*angle

    elif isinstance(q, (int, float)):
        return Quaternion(_log1p(q), 0, 0, 0)

    return NotImplemented

def pow(q: Quaternion or int or float, x: int or float) -> Quaternion:
    """Return q**x (q to the power of x), where x is a real exponent."""
    if (isinstance(q, Quaternion) and (isinstance(x, (int, float)))):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_pow(a, x), 0, 0, 0)

        theta = q.angle
        return (_pow(abs(q), x) * exp(q.unit_vector()*x*theta))

    elif isinstance(q, (int, float)):
        return Quaternion(_pow(q, x), 0, 0, 0)

    return NotImplemented

def rotate3d(point: tuple or list, angle: int or float,
             axis: tuple or list=(0, 0, 1), rounding: int=-1) -> tuple:
    """
    Takes a point in 3d space represented as a tuple or list of three
    (3) values and rotates it by an angle around a given axis vector.
    The axis of rotation is the z-axis [0, 0, 1] by default. For the
    point and axis parameters, if only one value is given, the value
    will be assumed to be an x-coordinate with the y- and z-coordinates
    equal to zero (0). If two values are given, they will be assumed to
    be x- and y-coordinates with the z-coordinate equal to zero (0).
    By default, this function does not round the result. If you wish to
    round the result to a number of decimal places, change the rounding
    argument to the number of decimal places you wish to round to.
    """
    if len(point) <= 3 and len(axis) <= 3:
        i, j, k = _makeListLen3(point)
        p = Quaternion(0, i, j, k)

        i2, j2, k2 = _makeListLen3(axis)
        u = Quaternion(0, i2, j2, k2)
        if abs(u) == 0:
            raise ValueError(
                'The axis to rotate around must be a nonzero vector.')
        if abs(u) != 1.0:
            u = u.versor  # Ensures u is a unit vector.
        q = _cos(angle/2) + u*_sin(angle/2)
        if abs(q) != 1.0:
            q = q.versor  # Same as above with u.
        p_prime = q * p * q.inverse()
        i_prime, j_prime, k_prime = p_prime.get_vector_components()
        if rounding != -1:
            i_prime, j_prime, k_prime = (
                round(i_prime, rounding), round(j_prime, rounding),
                round(k_prime, rounding))

        return (i_prime, j_prime, k_prime)

def _makeListLen3(i: Iterable[int or float]) -> list:
    """Makes sure points and axes of rotation have 3 coordinates."""
    if not isinstance(i, list):
        i = list(i)
    while len(i) < 3:
        i.append(0.0)
    return i
