# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 21:23:28 2021

@author: Zach Chartrand <zachartrand999@gmail.com>

Similar to the built-in module cmath, this module has definitions of
mathematical functions expanded to work with quaternions.
"""

__all__ = ['exp', 'log', 'log10', 'sqrt', 'rotate3d', 'pi', 'tau', 'e', 'inf',
           'infi', 'infj', 'infk', 'nan', 'nani', 'nanj', 'nank']

from math import (
    exp as _exp,
    cos as _cos,
    sin as _sin,
    log as _log,
    pi, tau, e, inf, nan,
)

from typing import Iterable

from quaternions import Quaternion

nani = Quaternion(0, float('nan'), 0, 0)
nanj = Quaternion(0, 0, float('nan'), 0)
nank = Quaternion(0, 0, 0, float('nan'))

infi = Quaternion(0, float('inf'), 0, 0)
infj = Quaternion(0, 0, float('inf'), 0)
infk = Quaternion(0, 0, 0, float('inf'))


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


def log(q: Quaternion or int or float, base: int or float=0) -> Quaternion:
    """Return the natural logarithm of a quaternion."""
    if isinstance(q, Quaternion):
        a = q.scalar
        if q.is_scalar():
            return Quaternion(_log(a), 0, 0, 0)

        angle = q.angle
        answer = _log(abs(q)) + q.unit_vector()*angle
        if base:
            answer = answer / _log(base)
        return answer

    elif isinstance(q, (int, float)):
        return Quaternion(_log(q), 0, 0, 0)

    return NotImplemented


def log10(q: Quaternion or int or float) -> Quaternion:
    """Return the base-10 logarithm of the quaternion."""
    return (log(q) / _log(10))


def sqrt(q: Quaternion or int or float) -> Quaternion:
    """Return the square root of q."""
    if q.is_scalar() and q.real < 0:
        raise ValueError("Negative real quaternions have an infinite number "
            + f"of square roots.\nThe square root of {q.real} is the sphere "
            + f"of radius {abs(q)**0.5:.4f}... centered at the origin.")

    return q.__pow__(0.5)


def _makeListLen3(i: Iterable[int or float]) -> list:
    """Makes sure points and axes of rotation have 3 coordinates."""
    if not isinstance(i, list):
        i = list(i)
    while len(i) < 3:
        i.append(0.0)
    return i


def rotate3d(point: tuple or list, angle: int or float,
             axis: tuple or list=(0, 0, 1), rounding: int=-1,
             degrees: bool=True) -> tuple:
    """
    Takes a point in 3d space represented as a tuple or list of three
    (3) values and rotates it by an angle around a given axis vector.
    The axis of rotation is the z-axis (0, 0, 1) by default. For the
    point and axis parameters, if only one value is given, the value
    will be assumed to be an x-coordinate with the y- and z-coordinates
    equal to zero (0). If two values are given, they will be assumed to
    be x- and y-coordinates with the z-coordinate equal to zero (0).
    By default, this function does not round the result. If you wish to
    round the result to a number of decimal places, change the rounding
    argument to the number of decimal places you wish to round to.
    The angle is set to be input in degrees. If you wish to use radians, set
    'degrees' equal to False.
    """
    if degrees:
        angle = angle % 360
        angle = angle * pi / 180
    else:
        angle = angle % tau

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
