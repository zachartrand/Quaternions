# -*- coding: utf-8 -*-
"""
Similar to the built-in module :py:mod:`cmath`, this module has
definitions of mathematical functions expanded to work with quaternions.
"""
# Metadata
# --------
# Created on Wed Apr 28 21:23:28 2021
# @author: Zach Chartrand <zachartrand999@gmail.com>

__all__ = [
    'isinf', 'isfinite', 'isnan', 'isclose', 'exp', 'log', 'log10', 
    'sqrt', 'cosh', 'cos', 'sinh', 'sin', 'tanh', 'tan', 'rotate3d', 
    'dot_product', 'cross_product', 'rotate_Euler', 'pi', 'tau', 'e', 
    'inf', 'infi', 'infj', 'infk', 'nan', 'nani', 'nanj', 'nank',
]

import sys
import math as _math
from math import pi, tau, e, inf, nan
import cmath as _cmath
from typing import Iterable, Tuple

from quaternions import Quaternion
from quaternions._misc import makeListLen3 as _makeListLen3

nani: Quaternion = Quaternion(0, nan, 0, 0)
nanj: Quaternion = Quaternion(0, 0, nan, 0)
nank: Quaternion = Quaternion(0, 0, 0, nan)

infi: Quaternion = Quaternion(0, inf, 0, 0)
infj: Quaternion = Quaternion(0, 0, inf, 0)
infk: Quaternion = Quaternion(0, 0, 0, inf)

_INV_LN10 = 0.4342944819032518  # Reciprocal of the natural log of 10.
_LN2 = 0.69314718055994530942  # Natural log of 2.
_HALF = 0.5
_TWO = 2.0
_SMALL_FLOAT = sys.float_info.min
_LARGE_FLOAT = 0.25*sys.float_info.max
_LOG_LARGE_FLOAT = _math.log(_LARGE_FLOAT)
_MANT_DIG = sys.float_info.mant_dig
_SCALE_DOWN_EXP = -_MANT_DIG//2 + 1
_SCALE_UP_EXP = 2*(-_SCALE_DOWN_EXP)


def isinf(q: Quaternion) -> bool:
    """
    Return ``True`` if any component of q is an infinity, 
    and ``False`` otherwise.
    """
    if isinstance(q, (int, float)):
        return _math.isinf(q)
    
    elif isinstance(q, Quaternion):
        return any([_math.isinf(component) for component in q])
    
    q_type = q.__class__.__qualname__
    raise TypeError(f"must be a Quaternion or real number, not {q_type}")


def isfinite(q: Quaternion) -> bool:
    """
    Return ``True`` if all components of q are finite, 
    and ``False`` otherwise.
    """
    if isinstance(q, (int, float)):
        return _math.isfinite(q)
    
    elif isinstance(q, Quaternion):
        return all([_math.isfinite(component) for component in q])
    
    q_type = q.__class__.__qualname__
    raise TypeError(f"must be a Quaternion or real number, not {q_type}")


def isnan(q: Quaternion) -> bool:
    """
    Return ``True`` if any component of q is a NaN, 
    and ``False`` otherwise.
    """
    if isinstance(q, (int, float)):
        return _math.isnan(q)
    
    elif isinstance(q, Quaternion):
        return any([_math.isnan(component) for component in q])

    q_type = q.__class__.__qualname__
    raise TypeError(f"must be a Quaternion or real number, not {q_type}")


def isclose(a, b, *, rel_tol=1e-09, abs_tol=1e-09) -> bool:
    """
    Determine whether two Quaternions are close in value.
    
    For the values to be considered close, the difference between 
    them must be smaller than at least one of the tolerances.

    -inf, inf and NaN behave similarly to the IEEE 754 Standard. 
    That is, NaN is not close to anything, even itself. 
    inf and -inf are only close to themselves.

    Parameters:
        a (Quaternion): The first Quaternion.
        b (Quaternion): The second Quaternion.
        rel_tol (float): maximum difference for being considered 
            "close", relative to the magnitude of the input values
        abs_tol (float): maximum difference for being considered 
            "close", regardless of the magnitude of the input values
    
    Returns:
        bool: ``True`` if a is close in value to b, 
            and ``False`` otherwise.
    """

    # Catch bad types.
    if not isinstance(a, (int, float, Quaternion)):
        a_type = a.__class__.__qualname__
        raise TypeError("must be a Quaternion or real number, "
                        f"not {a_type}")

    if not isinstance(b, (int, float, Quaternion)):
        b_type = b.__class__.__qualname__
        raise TypeError("must be a Quaternion or real number, "
                        f"not {b_type}")
    
    if not isinstance(rel_tol, (int, float)):
        rel_tol_type = rel_tol.__class__.__qualname__
        raise TypeError(f"must be real number, not {rel_tol_type}")

    if not isinstance(abs_tol, (int, float)):
        abs_tol_type = abs_tol.__class__.__qualname__
        raise TypeError(f"must be real number, not {abs_tol_type}")
    
    # Catch bad values for rel_tol and abs_tol.
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("tolerances must be non-negative")
    
    # Convert a and b to Quaternions if they are ints or floats.
    if isinstance(a, (int, float)):
        a = Quaternion(a, 0, 0, 0)
    
    if isinstance(b, (int, float)):
        b = Quaternion(b, 0, 0, 0)
    
    # Exact equality returns True.
    if a == b:
        return True
    
    # Infinities of the same sign are caught by the above equality.
    # This one catches infinities of differing signs and returns False.
    if (any([_math.isinf(component) for component in a]) 
        or any([_math.isinf(component) for component in b])):
        return False
    
    diff = (a - b).norm

    return (((diff <= rel_tol * b.norm) 
            or (diff <= rel_tol*a.norm)) 
            or diff <= abs_tol)


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
        if q.is_scalar():
            return Quaternion(_math.exp(q.real), q.i, q.j, q.k)  # Preserve signed zeroes.

        elif isfinite(q):
            vector = q.vector
            v_norm = q.vector_norm
            
            return (_math.exp(q.real)/v_norm)*(vector.cos_norm()*v_norm + vector*vector.sin_norm())
            
        else:  # Handle infinities and NaNs.
            try:
                imag = q.vector_norm
            
            except OverflowError:
                imag = _LARGE_FLOAT
            
            z = _cmath.exp(complex(q.real, imag))

            return q.from_complex(z)

    elif isinstance(q, (int, float)):
        return Quaternion(_math.exp(q), 0, 0, 0)

    # Raise a TypeError for complex numbers. Complex numbers should
    # either use the cmath function or be converted to a quaternion
    # before using this function.
    elif isinstance(q, complex):
        raise TypeError("Use the cmath module's exp() function or convert "
                        "the complex number to a quaternion.")

    # Raise a TypeError if the input is not a real number or quaternion.
    q_type = q.__class__.__qualname__
    raise TypeError("Cannot take the quaternion exponential of a "
                    f"{q_type}.")


def log(q: Quaternion or float, base: Quaternion or float = e) -> Quaternion:
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
    # ==================================================================
    if isinstance(q, (int, float)):
        return Quaternion(_math.log(q), 0.0, 0.0, 0.0)
    
    elif isinstance(q, Quaternion):
        # If the quaternion is real, use the normal logarithm function.
        if q.is_scalar():
            return Quaternion(_math.log(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        elif not isfinite(q):
            try:
                z = _cmath.log(complex(q.real, q.vector_norm))
            
            except OverflowError:
                q = _HALF*q
                z = _cmath.log(complex(q.real, q.vector_norm)) + _LN2
            
            result = q.from_complex(z)
            if base != e:
                result = result / log(base)
            
            return result
        
        real = q.log_norm()
        vector = (q.angle/q.vector_norm)*q.vector
        result = real + vector
        if base != e:
            result = result / log(base)
        
        return result

    q_type = q.__class__.__qualname__
    raise TypeError("Cannot take the quaternion logarithm of a "
                    f"{q_type}.")


def log10(q: Quaternion or float) -> Quaternion:
    """Return the base-10 logarithm of the quaternion."""
    return (log(q)*_INV_LN10)


def sqrt(q: Quaternion or float) -> Quaternion:
    """Return the square root of the quaternion."""
    sqrt_real = 0.0
    if isinstance(q, (int, float)):
        if q < 0.0:
            # Negative real numbers have an infinite number of 
            # square roots in the quaternion number system.
            real = q
            sqrt_real = abs(q)**_HALF
        elif q == 0.0:
            return Quaternion(0, 0, 0, 0)
        else:
            return Quaternion(_math.sqrt(q), 0, 0, 0)
    
    elif isinstance(q, Quaternion):
        if q.is_scalar():
            if q.real < 0:
                # Negative real numbers have an infinite number of 
                # square roots in the quaternion number system.
                real = q.real
                sqrt_real = (abs(q.real))**_HALF
            elif q.is_zero():
                return Quaternion(0, q.i, q.j, q.k)  # Preserve signed zeros.
            
            else:
                return Quaternion(_math.sqrt(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        else:  # Use cmath.sqrt function.
            # Find largest magnitude component and scale accordingly.
            max_component = max(q.abs_components())
            if max_component <= 4.0*_SMALL_FLOAT:
                q = _math.ldexp(1.0, _SCALE_UP_EXP)*q
                scale = _math.ldexp(1.0, _SCALE_DOWN_EXP)
            
            else:
                q = 0.25*q
                scale = _TWO
            
            imag = q.vector_norm
            z = _cmath.sqrt(complex(q.real, imag))
            
            return scale*q.from_complex(z)
    
    if sqrt_real:
        sqrt_string = str(sqrt_real)
        if len(str(sqrt_real)) > 8:
            sqrt_string = f"{sqrt_real:.6f}..."
        
        raise ValueError(
            "Negative real quaternions have an infinite number of square "
            + f"roots.\nThe square root of {real} is the sphere of "
            + f"vector quaternions of radius {sqrt_string} "
            + "centered at the origin.")
    
    else:
        q_type = q.__class__.__qualname__
        raise TypeError(
            f"Cannot take the Quaternion square root of a {q_type}")


def cosh(q: Quaternion) -> Quaternion:
    """Return the hyperbolic cosine of q."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.cosh(q), 0, 0, 0)

    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.cosh(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            imag = _math.hypot(*[_HALF*component for component in q.get_vector_components()])
            imag = _TWO*_math.fmod(imag, tau)
        
        z = _cmath.cosh(complex(q.real, imag))

        return q.from_complex(z)
    
    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion hyperbolic cosine of a {q_type}")


def cos(q: Quaternion or float) -> Quaternion:
    """Return the cosine of the quaternion."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.cos(q), 0, 0, 0)
    
    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.cos(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            raise OverflowError("math range error, vector norm too large")
        
        z = _cmath.cos(complex(q.real, imag))

        return q.from_complex(z)

    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion cosine of a {q_type}")


def sinh(q: Quaternion) -> Quaternion:
    """Return the hyperbolic sine of q."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.sinh(q), 0, 0, 0)

    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.sinh(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            imag = _math.hypot(*[_HALF*component for component in q.get_vector_components()])
            imag = _TWO*_math.fmod(imag, tau)
        
        z = _cmath.sinh(complex(q.real, imag))

        return q.from_complex(z)
    
    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion hyperbolic sine of a {q_type}")


def sin(q: Quaternion or float) -> Quaternion:
    """Return the sine of the quaternion."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.sin(q), 0, 0, 0)
    
    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.sin(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            raise OverflowError("math range error, vector norm too large")
        
        z = _cmath.sin(complex(q.real, imag))

        return q.from_complex(z)

    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion sine of a {q_type}")


def tanh(q: Quaternion) -> Quaternion:
    """Return the hyperbolic tangent of q."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.tanh(q), 0, 0, 0)

    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.tanh(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            imag = _math.hypot(*[_HALF*component for component in q.get_vector_components()])
            imag = _TWO*_math.fmod(imag, tau)
        
        z = _cmath.tanh(complex(q.real, imag))

        return q.from_complex(z)
    
    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion hyperbolic tangent of a {q_type}")


def tan(q: Quaternion or float) -> Quaternion:
    """Return the tangent of the quaternion."""
    if isinstance(q, (int, float)):
        return Quaternion(_math.tan(q), 0, 0, 0)
    
    elif isinstance(q, Quaternion):
        if q.is_scalar():
            return Quaternion(_math.tan(q.real), q.i, q.j, q.k)  # Preserve signed zeros.
        
        try:
            imag = q.vector_norm
        
        except OverflowError:
            # cmath.tan does not overflow with a very large imaginary part. 
            # In fact, it converges. _LOG_LARGE_FLOAT is large enough 
            # to reach that value.
            imag = _LOG_LARGE_FLOAT
        
        z = _cmath.tan(complex(q.real, imag))

        return q.from_complex(z)

    q_type = q.__class__.__qualname__
    raise TypeError(f"Cannot take the Quaternion tangent of a {q_type}")


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
        if _math.hypot(axis_i, axis_j, axis_k) < 1e-12:
            raise ValueError(
                'The axis of rotation must be a nonzero vector.')

        q = Quaternion.from_angle(
                angle*_HALF, (axis_i, axis_j, axis_k), degrees=degrees)
        
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
            angle parameters as degrees. Set this parameter to
            ``False`` to use angles in radians. Default is ``True``.
    """
    # yaw: rotate around z-axis
    # pitch: rotate around y'-axis
    # roll: rotate around x"-axis
    # TODO: Code all standards of Euler angles into this function.
    
    q_yaw = Quaternion.from_angle(yaw*_HALF, z_axis, degrees=degrees)
    y_axis = cross_product(z_axis, x_axis)
    q_pitch = Quaternion.from_angle(pitch*_HALF, y_axis, degrees=degrees)
    q_roll = Quaternion.from_angle(roll*_HALF, x_axis, degrees=degrees)
    q_total = q_yaw * q_pitch * q_roll

    q_point = Quaternion(0, *point)
    q_p_prime = q_total * q_point * q_total.inverse()

    return q_p_prime.get_vector_components()
