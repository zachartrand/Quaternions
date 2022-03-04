# -*- coding: utf-8 -*-
"""
Similar to the built-in module :py:mod:`cmath`, this module has
definitions of mathematical functions expanded to work with quaternions.
"""
# Metadata
# --------
# Created on Wed Apr 28 21:23:28 2021
# @author: Zach Chartrand <zachartrand999@gmail.com>

__all__ = ['exp', 'log', 'log10', 'sqrt', 'rotate3d', 'rotate_Euler',
           'dot_product', 'cross_product', 'pi', 'tau', 'e', 'inf',
           'infi', 'infj', 'infk', 'nan', 'nani', 'nanj', 'nank']

from math import (
    exp as _exp,
    cos as _cos,
    sin as _sin,
    log as _log,
    log1p as _log1p,
    hypot as _hypot,
    pi, tau, e, inf, nan,
)
from typing import Iterable, Tuple

from . import Quaternion
from .quaternions import _makeListLen3

nani: Quaternion = Quaternion(0, float('nan'), 0, 0)
nanj: Quaternion = Quaternion(0, 0, float('nan'), 0)
nank: Quaternion = Quaternion(0, 0, 0, float('nan'))

infi: Quaternion = Quaternion(0, float('inf'), 0, 0)
infj: Quaternion = Quaternion(0, 0, float('inf'), 0)
infk: Quaternion = Quaternion(0, 0, 0, float('inf'))


def exp(q: Quaternion or float) -> Quaternion:
    """Return the exponential of a quaternion."""
    # ==================================================================
    # The exponential (exp) of a quaternion (q) with a nonzero vector
    # component is given by the formula
    #
    #     exp(q) = exp(||q||) * (cos(||v||) + v*sin(||v||)/||v||)
    #
    # where ||q|| is the magnitude of the quaternion, v is the vector
    # part of the quaternion, and ||v|| is the magnitude of v. This
    # algorithm avoids this formula if the vector part is zero (0) or
    # if the input is an int or float.
    # ==================================================================
    if isinstance(q, Quaternion):
        a = q.real
        # If there is only a real component, skip all of the
        # complicated mathematics.
        if q.is_scalar():
            return Quaternion(_exp(a), 0, 0, 0)

        theta = q.vector_norm
        unit_vector = q.unit_vector()
        return _exp(a)*(_cos(theta) + unit_vector*_sin(theta))

    elif isinstance(q, (int, float)):
        return Quaternion(_exp(q), 0, 0, 0)

    # Raise a TypeError for complex numbers. Complex numbers should
    # either use the cmath function or be converted to a quaternion
    # before using this function.
    elif isinstance(q, complex):
        raise TypeError("Use the cmath module's exp() function or convert "
                      + "the complex number to a quaternion.")

    # Raise a TypeError if the input is not a real number or quaternion.
    raise TypeError("Cannot take the quaternion exponential of a "
                 + f"{q.__class__.__qualname__}.")


def log(q: Quaternion or float, base: float = e) -> Quaternion:
    """
    Return the logarithm of a quaternion to the given base.

    If the base is not specified, returns the natural logarithm
    (base *e*) of the quaternion.
    """
    # ==================================================================
    # The natural logarithm of a quaternion (q) is given by the formula
    #
    #     log(q) = log(||q||) + v*phi/||v||
    #
    # where ||q|| is the norm of the quaternion, v is the vector part
    # of the quaternion, ||v|| is the norm of the vector part of the
    # quaternion, and phi is the angle of the quaternion given by
    #
    #     phi = acos(a/||q||),
    #
    # where acos is the inverse cosine function and a is the scalar
    # part of the quaternion.
    #
    # Because log1p() is more accurate than log() for inputs close to
    # 1.0, when the norm of the quaternion is between 0.71 and 1.73
    # inclusive, log1p is automatically used. This range is used based
    # on the cmath module, which uses log1p for complex numbers with
    # magnitudes in this range.
    # ==================================================================
    if isinstance(q, Quaternion):
        real = q.real
        # If the quaternion is real, use the normal logarithm functions.
        if q.is_scalar():
            if 0.71 <= real and real <= 1.73:
                return Quaternion(_log1p(real-1), 0, 0, 0)
            else:
                return Quaternion(_log(real), 0, 0, 0)
        else:
            norm = q.norm
            if 0.71 <= norm and norm <= 1.73:
                # This takes advantage of the difference of perfect
                # squares to avoid potential loss of precision in
                # calculating the norm - 1 of the quaternion. This
                # formula is taken from the cmath log function,
                # modified for four (4) components.
                abs_components = [abs(real), abs(q.i), abs(q.j), abs(q.k)]
                max_component = max(abs_components)
                abs_components.remove(max_component)
                real = _log1p(
                    (max_component - 1) * (max_component + 1)
                    + abs_components[0]*abs_components[0]
                    + abs_components[1]*abs_components[1]
                    + abs_components[2]*abs_components[2]) / 2.0
            else:
                real = _log(norm)

            angle = q.angle
            answer = real + q.unit_vector()*angle
            if base != e:
                if 0.71 <= base and base <= 1.73:
                    answer = answer / _log1p(base-1)
                else:
                    answer = answer / _log(base)

            return answer

    elif isinstance(q, (int, float)):
        return log(Quaternion(q, 0.0, 0.0, 0.0))

    raise TypeError("Cannot take the quaternion logarithm of a "
                 + f"{q.__class__.__qualname__}.")


def log10(q: Quaternion or float) -> Quaternion:
    """Return the base-10 logarithm of the quaternion."""
    return (log(q) / _log(10))


def sqrt(q: Quaternion or float) -> Quaternion:
    """Return the square root of the quaternion."""
    # Negative real numbers have an infinite number of square roots in
    # the quaternion number system. These conditional statements check
    # if the input is a negative real number. If the input is not a
    # negative real number, this function returns the input raised to
    # 0.5. This should always return a Quaternion object.
    sqrt_real = 0.0
    if isinstance(q, Quaternion):
        if q.is_scalar() and q.real < 0:
            real = q.real
            sqrt_real = (-q.real)**0.5
    elif isinstance(q, (int, float)):
        if q < 0:
            real = q
            sqrt_real = (-q)**0.5
        else:
            q = Quaternion(q)

    if sqrt_real:
        if len(str(sqrt_real)) <= 4:
            sqrt_string = str(sqrt_real)
        else:
            sqrt_string = f"{sqrt_real:.6f}..."

        raise ValueError(
            "Negative real quaternions have an infinite number of square "
            + f"roots.\nThe square root of {real} is the sphere of radius "
            + f"{sqrt_string} centered at the origin.")

    return pow(q, 0.5)


def rotate3d(
        point: Iterable[float], angle: float,
        axis: Iterable[float] = (0.0, 0.0, 1.0), rounding: int = -1,
        degrees: bool = True) -> Tuple[float]:
    """
    Rotate a point around an axis.

    Take a point in 3d space represented as a tuple or list of three
    (3) values and rotate it by an angle around a given axis vector.

    Parameters:
        point: The point to rotate. The format for the coordinates is
            ``(x, y, z)``.
        angle: The angle of rotation. By default, angle is set to be
            input in degrees. See the **degrees** parameter if you want
            to use radians instead.
        axis: The axis to rotate the point around. By default, this is
            the z-axis ``(0, 0, 1)``.
        rounding: The number of decimal points the result will be
            rounded to. Default value is -1, which does not round the
            end result.
        degrees: When set to ``True``, this function interprets the
            parameter **angle** as degrees. Set this parameter to
            ``False`` to use angles in radians. Default is ``True``.

    For the point and axis parameters, if only one value is given, the
    value will be assumed to be an x-coordinate with the y- and
    z-coordinates equal to zero (0). If two values are given, they will
    be assumed to be x- and y-coordinates with the z-coordinate equal
    to zero (0).
    """
    if len(point) <= 3 and len(axis) <= 3:
        p = Quaternion(0, *_makeListLen3(point))

        axis_i, axis_j, axis_k = _makeListLen3(axis)
        if _hypot(axis_i, axis_j, axis_k) < 1e-12:
            raise ValueError(
                'The axis to rotate around must be a nonzero vector.')

        q = Quaternion.from_angle(
                angle*0.5, (axis_i, axis_j, axis_k), degrees=degrees)
        if abs(q) != 1.0:
            q = q.versor  # Ensures q is a unit vector.
        p_prime = q * p * q.inverse()
        if rounding >= 0:
            p_prime = round(p_prime, rounding)

        i_prime, j_prime, k_prime = p_prime.get_vector_components()

        return (i_prime, j_prime, k_prime)

    else:
        if len(point) > 3:
            problem = "'point'"
        elif len(axis) > 3:
            problem = "'axis'"
        else:
            problem = "'point' and 'axis'"
        raise IndexError(
            f"{problem} must be a tuple or list of three (3) values.")


def dot_product(
        vector1: Iterable[float], vector2: Iterable[float]) -> float:
    """
    Return the dot product of two vectors.

    Because this uses quaternions to calculate, this only works for
    vectors up to three (3) dimensions.
    """
    i1, j1, k1 = _makeListLen3(vector1)
    i2, j2, k2 = _makeListLen3(vector2)
    q1, q2 = Quaternion(0, i1, j1, k1), Quaternion(0, i2, j2, k2)

    return -(q1*q2).real


def cross_product(
    vector1: Iterable[float], vector2: Iterable[float]) -> Tuple[float]:
    """
    Return the cross product of two vectors.

    Because this uses quaternions to calculate, this only works for
    vectors up to three (3) dimensions.
    """
    i1, j1, k1 = _makeListLen3(vector1)
    i2, j2, k2 = _makeListLen3(vector2)
    q1, q2 = Quaternion(0, i1, j1, k1), Quaternion(0, i2, j2, k2)

    return (q1*q2).get_vector_components()

def rotate_Euler(
        point: Iterable[float],
        yaw: float, pitch: float, roll: float,
        x_axis: Iterable[float] = (1.0, 0.0, 0.0),
        z_axis: Iterable[float] = (0.0, 0.0, 1.0),
        degrees: bool = True) -> Tuple[float]:
    """
    Rotate a given point using Euler angles.

    This function uses the rotation convention of z-y'-x", rotating yaw,
    then pitch, then roll.

    Parameters:
        point: The point to rotate. The format for the coordinates is
            ``(x, y, z)``.
        yaw: The angle of rotation around the z-axis.
        pitch: The angle of rotation around the y'-axis. The y'-axis is
            the y-axis after the yaw rotation has been applied.
        roll: The angle of rotation around the x"-axis. The x"-axis is
            the x-axis after both the yaw and pitch rotations.
        x_axis: The initial x-axis of the coordinate system that
            **point** belongs to. Default value is ``(1, 0, 0)``.
        z_axis: The initial z-axis of the coordinate system that
            **point** belongs to. Default value is ``(0, 0, 1)``.
        degrees: When set to ``True``, this function interprets the
            parameter **angle** as degrees. Set this parameter to
            ``False`` to use angles in radians. Default is ``True``.
    """
    # yaw: rotate around z-axis
    # pitch: rotate around y'-axis
    # roll: rotate around x"-axis
    # TODO: Code all standards of Euler angles into this function.
    q_yaw = Quaternion.from_angle(yaw*0.5, z_axis)
    x_prime_axis = (q_yaw * Quaternion(0, *x_axis).versor
                    * q_yaw.inverse()).get_vector_components()
    y_prime_axis = cross_product(z_axis, x_prime_axis)
    q_pitch = Quaternion.from_angle(pitch*0.5, y_prime_axis)
    x_doubleprime_axis = (q_pitch * Quaternion(0, *x_prime_axis).versor
                    * q_pitch.inverse()).get_vector_components()
    q_roll = Quaternion.from_angle(roll*0.5, x_doubleprime_axis)

    q_point = Quaternion(0, *point)
    q_p_prime = (
        q_roll * (
            q_pitch * (
                q_yaw * q_point * q_yaw.inverse()
            ) * q_pitch.inverse()
        ) * q_roll.inverse()
    )

    return q_p_prime.get_vector_components()
