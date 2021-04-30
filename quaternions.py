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

from math import hypot as _hypot
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

            'a + bi + cj + dk',

        where a, b, c, and d are floats. If there is no real component,
        it is left out of the string, and has the format

            'bi + cj + dk'.
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
        if (self.real, self.i, self.j, self.k) == (0, 0, 0, 0):
            return False

        return True

    def __eq__(self, other):
        """Return self == other."""
        if isinstance(other, int) or isinstance(other, float):
            return (self.real, self.i, self.j, self.k) == (other, 0, 0, 0)
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

    def __neg__(self):
        return Quaternion(-self.real, -self.i, -self.j, -self.k)

    def __abs__(self):
        """Return abs(self)."""
        return _hypot(self.real, self.i, self.j, self.k)

    def __add__(self, other):
        """Return self+other."""
        if isinstance(other, int) or isinstance(other, float):
            real = self.real + other
            return Quaternion(real, self.i, self.j, self.k)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot add a Quaternion to a complex number. Make the '
                + 'complex number a Quaternion before adding.')
        elif isinstance(other, Quaternion):
            real = self.real + other.real
            i = self.i + other.i
            j = self.j + other.j
            k = self.k + other.k
            return Quaternion(real, i, j, k)
        else:
            raise TypeError('Cannot add a {} to a {}'.format(
                type(other).__name__, type(self).__name__))

    def __radd__(self, other):
        """Return other+self."""
        if isinstance(other, int) or isinstance(other, float):
            return self + other
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot add a complex number to a Quaternion. Make the '
                + 'complex number a Quaternion before adding.')

    def __sub__(self, other):
        """Return self-other."""
        if isinstance(other, int) or isinstance(other, float):
            real = self.real - other
            return Quaternion(real, self.i, self.j, self.k)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot subtract a complex number from a Quaternion. Make the '
                + 'complex number a Quaternion before subtracting.')
        elif isinstance(other, Quaternion):
            real = self.real - other.real
            i = self.i - other.i
            j = self.j - other.j
            k = self.k - other.k
            return Quaternion(real, i, j, k)
        else:
            raise TypeError('Cannot add a {} to a {}'.format(
                type(other).__name__, type(self).__name__))

    def __rsub__(self, other):
        """Return other-self."""
        if isinstance(other, int) or isinstance(other, float):
            real = other - self.real
            return Quaternion(real, -self.i, -self.j, -self.k)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot subtract a Quaternion from a complex number. Make the '
                + 'complex number a Quaternion before subtracting.')

    def __mul__(self, other):
        """
        Return self*other.

        If the Quaternion is multiplied by an int or float, the
        int/float is distributed to each component.

        A Quaternion multiplied by a complex number will raise a
        TypeError, as it is ambiguous which vector component the
        imaginary component should be treated as.
        """
        if isinstance(other, int) or isinstance(other, float):
            return Quaternion(
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
            raise TypeError('Cannot multiply a {} by a {}'.format(
                type(other).__name__, type(self).__name__))

        return Quaternion(real, i, j, k)

    def __rmul__(self, other):
        """Return other*self."""
        if isinstance(other, int) or isinstance(other, float):
            return self * other
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot multiply a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before multiplying.')

    def conjugate(self):
        """Returns the complex conjugate of self."""
        return Quaternion(self.real, -self.i, -self.j, -self.k)

    def inverse(self):
        """Return 1/self."""
        q_inverse = self.conjugate() / (
            self.real**2 + self.i**2 + self.j**2 + self.k**2)
        return q_inverse

    def __truediv__(self, other):
        """
        Return self/other.

        For division q1/q2, this assumes the order of multiplication
        is q1 * 1/q2.
        """
        if isinstance(other, int) or isinstance(other, float):
            return Quaternion(
                self.real/other, self.i/other, self.j/other, self.k/other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot divide a Quaternion by a complex number. Make the '
                + 'complex number a Quaternion before dividing.')
        elif isinstance(other, Quaternion):
            return self * other.inverse()

    def __rtruediv__(self, other):
        """Return other/self."""
        if isinstance(other, int) or isinstance(other, float):
            return self.inverse() * other
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot divide a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before dividing.')

    @property
    def scalar(self) -> float:
        """Return the real part of the quaternion."""
        return self.real

    @property
    def vector(self) -> List[float]:
        """Return the vector part of the quaternion."""
        return Quaternion(0, self.i, self.j, self.k)

    @property
    def unit_vector(self):
        """
        Return the vector part of the quaternion normalized to a
        magnitude of one (1).
        """
        if (self.i, self.j, self.k) == (0, 0, 0):
            return Quaternion(0, 0, 0, 0)
        else:
            v = Quaternion(0, self.i, self.j, self.k)
            return v/abs(v)

    def get_vector_components(self):
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
        if (self.i, self.j, self.k) == (0, 0, 0):
            return True
        return False

    def is_vector(self) -> bool:
        """
        Returns True if the scalar part is zero and at least one of
        the vector components is nonzero. Otherwise, returns False.
        """
        if self.real == 0 and (self.i != 0 or self.j != 0 or self.k != 0):
            return True
        return False

    def is_complex(self):
        """
        Returns True if only one of the i, j, and k components is
        nonzero. Otherwise, returns False.
        """
        if (self.i, self.j, self.k) != (0, 0, 0):
            if (0, 0) in (
                    (self.i, self.j), (self.j, self.k), (self.i, self.k)):
                return True

        return False

    def get_imag(self):
        """
        Returns the imaginary component of the quaternion if only one
        of the imaginary components is nonzero. Otherwise,
        returns None.
        """
        if self.is_complex():
            for component in (self.i, self.j, self.k):
                if component != 0:
                    return component
        elif self.is_scalar():
            return 0.0

        else:
            return None
