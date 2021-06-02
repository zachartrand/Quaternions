# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 15:36:03 2021

@author: Zach Chartrand <zachartrand999@gmail.com>

This module contains the Quaternion class, a class based on the
quaternion numbers. This class is designed to properly perform
quaternion arithmetic with other quaternions and int and float type
numbers. For complex numbers, this class does not assume which
quaternion vector the imaginary component belongs to, so in order to do
arithmetic with complex numbers, they must be converted into a
quaternion first.
"""

__all__ = ['Quaternion']

from math import hypot as _hypot, ceil as _ceil, atan2 as _atan2
from typing import List

class Quaternion():
    """
    Creates a quaternion. Quaternions are an expansion of the complex
    numbers, where there are four (4) components--the real component,
    also known as the scalar part, and the imaginary components, which
    together are known as the vector part. The vector part is made up
    of three (3) components whose unit values are i, j, and k. The
    rules for these values are as follows:

        i**2 == j**2 == k**2 == -1
        i*j == -j*i == k
        j*k == -k*j == i
        k*i == -i*k == j,

    which leads to the following statement:

        i*j*k == -1.
    """
    def __init__(self, real: int or float=0, i_component: int or float=0,
            j_component: int or float=0, k_component: int or float=0):
        self.real = float(real)
        self.i = float(i_component)
        self.j = float(j_component)
        self.k = float(k_component)

    def __repr__(self):
        """
        Return repr(self).

        Returns a string of the quaternion in the form

            '(a + bi + cj + dk)',

        where a, b, c, and d are floats. If there is no real component,
        it is left out of the string, and has the format

            '(bi + cj + dk)'.
        """
        i_str = str(self.i) + 'i'
        j_str = str(self.j) + 'j'
        k_str = str(self.k) + 'k'
        if self.real:
            q_str = ' + '.join([str(self.real), i_str, j_str, k_str])
        else:
            q_str = ' + '.join([i_str, j_str, k_str])

        return f'({q_str})'

    def __bool__(self):
        """self != 0"""
        if (self.real, self.i, self.j, self.k) == (0.0, 0.0, 0.0, 0.0):
            return False

        return True

    def __eq__(self, other):
        """Return self == other."""
        if isinstance(other, (int, float)):
            return (self.real, self.i, self.j, self.k) == (other, 0.0, 0.0, 0.0)
        elif isinstance(other, complex):
            return (self.real, self.get_imag()) == (
                other.real, other.imag)
        elif isinstance(other, Quaternion):
            return (self.real, self.i, self.j, self.k) == (
                other.real, other.i, other.j, other.k)

        return False

    def __hash__(self):
        """
        Return hash(self).

        Works so that scalar quaternions have the same hash values as
        int and float. Quaternions that have two of their imaginary
        terms equal to zero (0) have the same hash values as complex
        numbers.
        """
        if self.is_scalar():
            return hash(self.real)
        elif self.is_complex():
            return hash(complex(self.real, self.get_imag()))
        else:
            return hash((self.real, self.i, self.j, self.k))

    def __iter__(self):
        """Returns an iterator of the components of the quaternion."""
        return iter([self.real, self.i, self.j, self.k])

    def __reversed__(self):
        """
        Returns an iterator of the components of the quaternion
        in reverse order.
        """
        return reversed([self.real, self.i, self.j, self.k])

    def __pos__(self):
        """Return +self."""
        return self

    def __neg__(self):
        """Return -self."""
        return Quaternion(-self.real, -self.i, -self.j, -self.k)

    def __abs__(self):
        """
        Return abs(self).

        Returns the magnitude of the quaternion.
        """
        return _hypot(self.real, self.i, self.j, self.k)

    def __int__(self):
        """
        Return int(self).

        If the quaternion is scalar, this returns the real component as
        an integer.
        """
        if self.is_scalar():
            return int(self.real)

        return NotImplemented

    def __float__(self):
        """
        Return float(self).

        If the quaternion is scalar, returns the scalar component.
        """
        if self.is_scalar():
            return self.real

        return NotImplemented

    def __round__(self, ndigits=None):
        """
        Return round(self, ndigits).

        Rounds each component of the quaternion to the nearest integer.
        Each component is returned as a float.
        """
        return Quaternion(
            round(self.real, ndigits), round(self.i, ndigits),
            round(self.j, ndigits), round(self.k, ndigits))

    def __floor__(self):
        """
        Return math.floor(self)

        Returns the quaternion with all of its components rounded down
        to the nearest integer.
        """
        return Quaternion(
            self.real // 1, self.i // 1, self.j // 1, self.k // 1)

    def __ceil__(self):
        """
        Return math.ceil(self).

        Returns the quaternion with all of its components rounded up
        to the nearest integer.
        """
        return Quaternion(
            _ceil(self.real), _ceil(self.i),
            _ceil(self.j), _ceil(self.k))

    def __add__(self, other):
        """
        Return self+other.

        If other is an int or float, it gets added to the scalar part
        of the quaternion. If other is a quaternion, the quaternions
        are added component-by-component.
        """
        if isinstance(other, (int, float)):
            real = self.real + other
            i, j, k = self.i, self.j, self.k
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot add a Quaternion to a complex number. Make the '
                + 'complex number a Quaternion before adding.')
        elif isinstance(other, Quaternion):
            real = self.real + other.real
            i = self.i + other.i
            j = self.j + other.j
            k = self.k + other.k
        else:
            return NotImplemented

        return Quaternion(real, i, j, k)

    def __radd__(self, other):
        """Return other+self."""
        if isinstance(other, (int, float)):
            return self.__add__(other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot add a complex number to a Quaternion. Make the '
                + 'complex number a Quaternion before adding.')

        return NotImplemented

    def __sub__(self, other):
        """
        Return self-other.

        If other is an int or float, it is subtracted from the scalar
        part of the quaternion. If other is a quaternion, the second
        quaternion is subtracted by the first component by component.
        """
        if isinstance(other, (int, float)):
            real = self.real - other
            i, j, k = self.i, self.j, self.k
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot subtract a complex number from a Quaternion. Make the '
                + 'complex number a Quaternion before subtracting.')
        elif isinstance(other, Quaternion):
            real = self.real - other.real
            i = self.i - other.i
            j = self.j - other.j
            k = self.k - other.k
        else:
            return NotImplemented

        return Quaternion(real, i, j, k)

    def __rsub__(self, other):
        """
        Return other-self.

        If other is an int or float, the scalar part of the quaternion
        is subtracted from other, and the signs of the vector
        components are switched.
        """
        if isinstance(other, (int, float)):
            real = other - self.real
            return Quaternion(real, -self.i, -self.j, -self.k)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot subtract a Quaternion from a complex number. Make the '
                + 'complex number a Quaternion before subtracting.')

        return NotImplemented

    def __mul__(self, other):
        """
        Return self*other.

        If the Quaternion is multiplied by an int or float, the
        int/float is distributed to each component. If other is a
        quaternion, this multiplies following the rules of quaternion
        multiplication. Multiplication is not commutative; q1 * q2
        usually does not have the same result as q2 * q1.
        """
        if isinstance(other, (int, float)):
            real, i, j, k = (
                other*self.real, other*self.i, other*self.j, other*self.k)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot multiply a Quaternion by a complex number. Make the '
                + 'complex number a Quaternion before multiplying.')
        elif isinstance(other, Quaternion):
            real = (self.real*other.real - self.i*other.i - self.j*other.j
                    - self.k*other.k)
            i = (self.real*other.i + self.i*other.real + self.j*other.k
                 - self.k*other.j)
            j = (self.real*other.j + self.j*other.real + self.k*other.i
                 - self.i*other.k)
            k = (self.real*other.k + self.k*other.real + self.i*other.j
                 - self.j*other.i)
        else:
            return NotImplemented

        return Quaternion(real, i, j, k)

    def __rmul__(self, other):
        """Return other*self."""
        if isinstance(other, (int, float)):
            return self.__mul__(other)  # Scalar multiplication commutes.
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot multiply a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before multiplying.')

        return NotImplemented

    def conjugate(self):
        """
        Returns the conjugate of self. This is analogous to the
        complex conjugate, reversing the signs of the vector
        components.
        """
        return Quaternion(self.real, -self.i, -self.j, -self.k)

    def inverse(self):
        """
        Return 1/self.

        Returns the inverse of the quaternion. The inverse of a
        quaternion is defined as the conjugate divided by the norm
        squared:
            q.inverse() = q.conjugate()/(q.norm)**2.
        """
        # =====================================================================
        # Added an algorithm similar to Python's algorithm for complex numbers
        # to avoid overflow and underflow errors. Using this algorithm avoids
        # squaring numbers and scales the components by the largest of the
        # components, making the base numbers to work with between 0 and 1
        # in magnitude.
        # =====================================================================
        max_component = max(abs(self.real), abs(self.i), abs(self.j), abs(self.k))
        real_ratio = self.real / max_component
        i_ratio = self.i / max_component
        j_ratio = self.j / max_component
        k_ratio = self.k / max_component
        denom = (real_ratio * self.real + i_ratio * self.i
                 + j_ratio * self.j + k_ratio * self.k)

        q_inverse = (
            Quaternion(real_ratio, -i_ratio, -j_ratio, -k_ratio)
            / denom)
        return q_inverse

    def __floordiv__(self, other):
        """
        Return self // other.

        If other is an int or float, returns the quaternion with each
        component floor divided by other.
        """
        if isinstance(other, (int, float)):
            return Quaternion(
                self.real // other, self.i // other, self.j // other,
                self.k // other)
        else:
            return NotImplemented

    def __truediv__(self, other):
        """
        Return self/other.

        For division q1/q2, this assumes the order of multiplication
        is q1 * 1/q2. To left-multiply the denominator, enter 1/q2 * q1.
        """
        if isinstance(other, (int, float)):
            return Quaternion(
                self.real/other, self.i/other, self.j/other, self.k/other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot divide a Quaternion by a complex number. Make the '
                + 'complex number a Quaternion before dividing.')
        elif isinstance(other, Quaternion):
            return self.__mul__(other.inverse())
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        """Return other/self."""
        if isinstance(other, (int, float)):
            return self.inverse().__mul__(other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot divide a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before dividing.')

        return NotImplemented

    def __complex__(self):
        """
        Return complex(self).

        If only one of the vector components is nonzero, this returns the
        quaternion as a complex number.
        """
        if self.is_complex():
            real = self.real
            imag = self.get_imag()
            return complex(real, imag)
        return NotImplemented

    @property
    def scalar(self) -> float:
        """Return the real part of the quaternion."""
        return self.real

    @property
    def vector(self):
        """Return the vector part of the quaternion."""
        return Quaternion(0, self.i, self.j, self.k)

    @property
    def norm(self):
        """Returns the norm (magnitude) of the quaternion."""
        return abs(self)

    @property
    def vector_norm(self):
        """Returns the norm of the vector part of the quaternion."""
        return abs(self.vector)

    @property
    def angle(self):
        """
        Returns the angle of the quaternion.

        Quaternions can be expressed as

            norm*(cos(theta) + u*sin(theta))

        where u is a 3D unit vector. This returns the angle theta from
        this expression.
        """
        return _atan2(abs(self.vector), self.real)

    @property
    def versor(self):
        """Returns the quaternion normalized to a magnitude of one (1)."""
        versor = self
        while abs(versor) != 1.0:  # Doesn't always divide to norm 1 on the
                                   # first division.
            versor = versor/abs(versor)
        return versor

    def unit_vector(self):
        """
        Return the vector part of the quaternion normalized to a
        magnitude of one (1). Returns the zero quaternion if the
        magnitude of the quaternion is zero (0).
        """
        if (self.i, self.j, self.k) == (0.0, 0.0, 0.0):
            return Quaternion(0, 0, 0, 0)
        else:
            v = Quaternion(0, self.i, self.j, self.k)
            return v.versor

    def unit_quaternion(self):
        """
        Return the quaternion normalized to unity (1).

        If the quaternion is a zero (0) quaternion, returns None.
        """
        if abs(self) != 0.0:
            return self.versor

        return None

    def get_vector_components(self) -> List[float]:
        """Return the vector components of the Quaternion as a list
        formatted as

            [i, j, k].
        """
        return [self.i, self.j, self.k]

    def is_scalar(self) -> bool:
        """
        Returns True if the vector components all equal zero.
        Otherwise, returns False.
        """
        if (self.i, self.j, self.k) == (0.0, 0.0, 0.0):
            return True

        return False

    def is_vector(self) -> bool:
        """
        Returns True if the scalar part is zero and at least one of
        the vector components is nonzero. Otherwise, returns False.
        """
        if self.real == 0.0 and (
                self.i != 0.0 or self.j != 0.0 or self.k != 0.0):
            return True

        return False

    def is_complex(self) -> bool:
        """
        Returns True if only one of the i, j, and k components is
        nonzero. Otherwise, returns False.
        """
        if (self.i, self.j, self.k) != (0.0, 0.0, 0.0):
            if (0.0, 0.0) in (
                    (self.i, self.j), (self.j, self.k), (self.i, self.k)):
                return True

        return False

    def get_imag(self) -> float:
        """
        Returns the imaginary component of the quaternion if only one
        of the imaginary components is nonzero. If the quaternion is
        scalar, this returns 0.0. Otherwise, returns None.
        """
        if self.is_complex():
            for component in (self.i, self.j, self.k):
                if component != 0.0:
                    return component
        elif self.is_scalar():
            return 0.0
        else:
            return None
