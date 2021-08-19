# -*- coding: utf-8 -*-
"""
This module contains the Quaternion class, a class based on the
quaternion numbers. This class is designed to properly perform
quaternion arithmetic with other quaternions and int and float type
numbers. For complex numbers, this class does not assume which
quaternion vector the imaginary component belongs to, so in order to do
arithmetic with complex numbers, they must be converted into a
quaternion first.
"""
# Metadata
# --------
# Created on Mon Apr 26 15:36:03 2021
# @author: Zach Chartrand <zachartrand999@gmail.com>

from __future__ import annotations

__all__ = ['Quaternion']

from math import (
    hypot as _hypot,
    floor as _floor,
    ceil as _ceil,
    atan2 as _atan2,
    copysign as _copysign,
    exp as _exp,
    cos as _cos,
    sin as _sin,
    pi as _pi,
    log as _log,
)

from typing import Tuple, Iterable, Iterator


def _makeListLen3(i: Iterable[float]) -> list:
    """Makes sure points and axes of rotation have 3 coordinates."""
    if not isinstance(i, list):
        i = list(i)

    if len(i) > 3:
        i = i[:3]
    elif len(i) < 3:
        while len(i) < 3:
            i.append(0.0)
    return i


class Quaternion():
    """
    Quaternions are an expansion of the complex numbers, where there are
    four (4) components--the real component, also known as the scalar part,
    and the imaginary components, which together are known as the vector part.
    The vector part is made up of three (3) components whose unit values
    are `i`, `j`, and `k`. The rules for these values are as follows:

      :math:`i^2 = j^2 = k^2 = -1`

      :math:`jk = -kj = i`

      :math:`ki = -ik = j`

      :math:`ij = -ji = k`,

    which leads to the following statement:

      :math:`ijk = -1`.

    The descriptions will reference a quaternion of the form
    :math:`a + bi + cj + dk`, where :math:`a`, :math:`b`, :math:`c`, and :math:`d`
    are real numbers.

    Parameters:
        real: The real component (:math:`a`) of the quaternion.
        i_component: The i component (:math:`b`) of the quaternion.
        j_component: The j component (:math:`c`) of the quaternion.
        k_component: The k component (:math:`d`) of the quaternion.

    Each component can be returned by calling the attribute of the same name.

    Example:
        >>> q = Quaternion(1, -2, -3, 4)
        >>> print(q)
        (1 - 2i - 3j + 4k)
        >>> q.real
        1.0
        >>> q.i
        -2.0
        >>> q.j
        -3.0
        >>> q.k
        4.0

    """
    ## Initialization ##
    ####################
    def __init__(self, real: float = 0.0, i_component: float = 0.0,
                 j_component: float = 0.0, k_component: float = 0.0) -> None:
        self._real: float = float(real)
        self._i: float = float(i_component)
        self._j: float = float(j_component)
        self._k: float = float(k_component)

    # Make all the components read-only attributes.
    @property
    def real(self) -> float:
        """The real component of the Quaternion."""
        return self._real

    @property
    def i(self) -> float:
        """The i component of the Quaternion."""
        return self._i

    @property
    def j(self) -> float:
        """The j component of the Quaternion."""
        return self._j

    @property
    def k(self) -> float:
        """The k component of the Quaternion."""
        return self._k

    ### Magic methods. ###
    ######################

    ## Object data magic methods. ##
    ################################
    def __complex__(self) -> complex:
        """
        Return complex(self).

        If only one of the vector components is nonzero, return the
        quaternion as a complex number.
        """
        if self.is_complex():
            real = self.real
            imag = self.get_imag()
            return complex(real, imag)

        elif self.is_scalar():
            real, imag = self.real, 0.0
            return complex(real, imag)

        raise ValueError(f"{self.__class__.__qualname__} is not scalar or complex.")

    def __float__(self) -> float:
        """
        Return float(self).

        If the quaternion is scalar, return the scalar component.
        """
        if self.is_scalar():
            return self.real

        raise ValueError(f"{self.__class__.__qualname__} is not scalar.")

    def __hash__(self) -> int:
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

    def __int__(self) -> int:
        """
        Return int(self).

        If the quaternion is scalar, return the real component as
        an integer.
        """
        if self.is_scalar():
            return int(self.real)

        raise ValueError(f"{self.__class__.__qualname__} is not scalar.")

    def __iter__(self) -> Iterator[float]:
        """
        Return iter(self).

        Return an iterator of the components of the quaternion.

        The components are in the order of

            real, i, j, k
        """
        return iter([self.real, self.i, self.j, self.k])

    def __reversed__(self) -> Iterator[float]:
        """
        Return reversed(self).

        Return an iterator of the components of the quaternion
        in reverse order:

            k, j, i, real
        """
        return reversed([self.real, self.i, self.j, self.k])

    def __repr__(self) -> str:
        """Return repr(self)."""
        return (f"{self.__class__.__qualname__}("
              + f"{self.real}, {self.i}, {self.j}, {self.k})")

    def __str__(self) -> str:
        """
        Return str(self).

        Return a string of the quaternion in the form

        ``'(a + bi + cj + dk)'``,

        where a, b, c, and d are floats. If there is no real component,
        it is left out of the string, and has the format

            ``'(bi + cj + dk)'``.
        """
        # Helper functions.
        def _sign(x: float) -> str:
            """Return '+' or '-' based on whether x is positive or negative."""
            if _copysign(1.0, x) == -1.0:
                return "-"
            else:
                return "+"

        def _num_to_str(x: float) -> str:
            """
            Return a string of x as an integer if x is a positive or
            negative whole number, otherwise return a float string.
            """
            if x.is_integer():
                return str(int(x))
            else:
                return str(x)

        j_str = "".join([_sign(self.j), " ", _num_to_str(self.j.__abs__()), 'j'])
        k_str = "".join([_sign(self.k), " ", _num_to_str(self.k.__abs__()), 'k'])
        if self.real:
            i_str = "".join([_sign(self.i), " ", _num_to_str(self.i.__abs__()), 'i'])
            q_str = ' '.join([_num_to_str(self.real), i_str, j_str, k_str])
        else:
            i_str = "".join([_num_to_str(self.i), 'i'])
            q_str = ' '.join([i_str, j_str, k_str])

        return f'({q_str})'

    ## Boolean magic methods ##
    ###########################
    def __bool__(self) -> bool:
        """self != 0"""
        if (self.real, self.i, self.j, self.k) == (0.0, 0.0, 0.0, 0.0):
            return False

        return True

    def __eq__(self, other: Quaternion or float or complex) -> bool:
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

    def __ge__(self, other):
        return NotImplemented

    def __gt__(self, other):
        return NotImplemented

    def __le__(self, other):
        return NotImplemented

    def __lt__(self, other):
        return NotImplemented

    ## Mathematical magic methods with one mandatory input. ##
    ##########################################################
    def __abs__(self) -> float:
        """
        Return abs(self).

        Return the magnitude of the quaternion.
        """
        return _hypot(self.real, self.i, self.j, self.k)

    def __ceil__(self) -> Quaternion:
        """
        Return math.ceil(self).

        Return the quaternion with all of its components rounded up
        to the nearest integer.
        """
        return Quaternion(
            _ceil(self.real), _ceil(self.i), _ceil(self.j), _ceil(self.k))

    def __floor__(self) -> Quaternion:
        """
        Return math.floor(self).

        Return the quaternion with all of its components rounded down
        to the nearest integer.
        """
        return Quaternion(
            _floor(self.real), _floor(self.i), _floor(self.j), _floor(self.k))

    def __pos__(self) -> Quaternion:
        """Return +self."""
        return self

    def __neg__(self) -> Quaternion:
        """Return -self."""
        return Quaternion(-self.real, -self.i, -self.j, -self.k)

    def __round__(self, ndigits: None or int = None) -> Quaternion:
        """
        Return round(self, ndigits).

        Rounds each component of the quaternion to the nearest integer.
        Each component is returned as a float.
        """
        return Quaternion(
            round(self.real, ndigits), round(self.i, ndigits),
            round(self.j, ndigits), round(self.k, ndigits))

    ## Mathematical magic methods with two mandatory inputs. ##
    ###########################################################
    def __add__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return self + other.

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

    def __radd__(self, other: Quaternion or float or complex) -> Quaternion:
        """Return other + self."""
        if isinstance(other, (int, float)):
            return self.__add__(other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot add a complex number to a Quaternion. Make the '
                + 'complex number a Quaternion before adding.')

        return NotImplemented

    def __sub__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return self - other.

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

    def __rsub__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return other - self.

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

    def __mul__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return self * other.

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

    def __rmul__(self, other: Quaternion or float or complex) -> Quaternion:
        """Return other * self."""
        if isinstance(other, (int, float)):
            return self.__mul__(other)  # Scalar multiplication commutes.
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot multiply a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before multiplying.')

        return NotImplemented

    def __floordiv__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return self // other.

        If other is an int or float, return the quaternion with each
        component floor divided by other.
        """
        if isinstance(other, (int, float)):
            return Quaternion(
                self.real // other, self.i // other, self.j // other,
                self.k // other)
        else:
            return NotImplemented

    def __truediv__(self, other: Quaternion or float or complex) -> Quaternion:
        """
        Return self/other.

        For division q1/q2, this assumes the order of multiplication
        is q1 * 1/q2. To left-multiply the denominator, enter 1/q2 * q1 or
        q2.inverse() * q1.
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

    def __rtruediv__(self, other: Quaternion or float or complex) -> Quaternion:
        """Return other/self."""
        if isinstance(other, (int, float)):
            return self.inverse().__mul__(other)
        elif isinstance(other, complex):
            raise TypeError(
                'Cannot divide a complex number by a Quaternion. Make the '
                + 'complex number a Quaternion before dividing.')

        return NotImplemented

    def __pow__(self, other: Quaternion or float,
                mod: None = None) -> Quaternion:
        """
        Return pow(self, other, mod).

        Return self to the power of other, where other is a real exponent.
        """
        if mod is not None:
            raise ValueError(
                f"{self.__class__.__qualname__}s are not set up "
              + "for modular arithmetic.")
        if (isinstance(other, int)
                or (isinstance(other, float) and other.is_integer())):
            x = int(other)
            if x == 1:
                return self
            elif x == -1:
                return self.inverse()
            q = Quaternion(1)
            if x > 1:
                for _ in range(x):
                    q = q.__mul__(self)
            elif x < -1:
                for _ in range(-x):
                    q = q.__mul__(self.inverse())

            return q

        elif isinstance(other, float):
            if self.is_scalar():
                return Quaternion(pow(self.real, other), 0, 0, 0)

            theta = self.angle
            return (
                pow(self.__abs__(), other) * (
                    _cos(other*theta)
                    + self.unit_vector().__mul__(_sin(other*theta))
                )
            )

        elif isinstance(other, Quaternion):
            if other.is_scalar():
                return self.__pow__(other.real)

            ln = _log(self.norm) + self.unit_vector()*self.angle
            q = ln * other
            theta = q.vector_norm % (_pi*2)
            pow_q = (
                _exp(q.real) * (_cos(theta) + q.unit_vector()*_sin(theta)))

            return pow_q

        return NotImplemented

    def __mod__(self, other):
        return NotImplemented

    ## Normal methods ##
    ####################
    def conjugate(self) -> Quaternion:
        """
        Return the conjugate of self. This is analogous to the
        complex conjugate, reversing the signs of the vector
        components.
        """
        return Quaternion(self.real, -self.i, -self.j, -self.k)

    def get_imag(self) -> float:
        """
        Return the imaginary component of the quaternion if only one
        of the imaginary components is nonzero. If the quaternion is
        scalar, return ``0.0``. Otherwise, return ``None``.
        """
        if self.is_complex():
            for component in (self.i, self.j, self.k):
                if component != 0.0:
                    return component
        elif self.is_scalar():
            return 0.0
        else:
            return None

    def get_vector_components(self) -> Tuple[float]:
        """
        Return the vector components of the Quaternion as a tuple
        formatted as ``(i, j, k)``.
        """
        return (self.i, self.j, self.k)

    def inverse(self) -> Quaternion:
        """
        Return 1/self.

        Return the inverse of the quaternion. The inverse of a
        quaternion is defined as the conjugate divided by the norm
        squared::

            q.inverse() = q.conjugate()/(q.norm)**2

        """
        # =====================================================================
        # Added an algorithm similar to Python's algorithm for complex numbers
        # to avoid overflow and underflow errors. Using this algorithm avoids
        # squaring numbers and scales the components by the largest of the
        # components, making the base numbers to work with between 0 and 1
        # in magnitude.
        # =====================================================================
        max_component = max(
            self.real.__abs__(), self.i.__abs__(), self.j.__abs__(),
            self.k.__abs__())
        real_ratio = self.real / max_component
        i_ratio = self.i / max_component
        j_ratio = self.j / max_component
        k_ratio = self.k / max_component
        denom = (real_ratio * self.real + i_ratio * self.i
                 + j_ratio * self.j + k_ratio * self.k)

        q_inverse = (
            Quaternion(real_ratio, -i_ratio, -j_ratio, -k_ratio) / denom)

        return q_inverse

    def unit_quaternion(self) -> Quaternion:
        """
        Return the quaternion normalized to magnitude one (1).

        If the quaternion is a zero (0) quaternion, return the zero quaternion.
        """
        if self.__abs__() != 0.0:
            return self.versor

        return Quaternion(0)

    def unit_vector(self) -> Quaternion:
        """
        Return the vector part of the quaternion normalized to a
        magnitude of one (1). Return the zero quaternion if the
        magnitude of the quaternion is zero (0).
        """
        if (self.i, self.j, self.k) == (0.0, 0.0, 0.0):
            return Quaternion(0, 0, 0, 0)
        else:
            v = Quaternion(0, self.i, self.j, self.k)
            return v.versor

    ## Boolean methods ##
    #####################
    def is_complex(self) -> bool:
        """
        Return ``True`` if only one of the *i*, *j*, and *k* components is
        nonzero. Otherwise, return ``False``.
        """
        if (self.i, self.j, self.k) != (0.0, 0.0, 0.0):
            if (0.0, 0.0) in (
                    (self.i, self.j), (self.j, self.k), (self.i, self.k)):
                return True

        return False

    def is_scalar(self) -> bool:
        """
        Return ``True`` if the vector components all equal zero.
        Otherwise, return ``False``.
        """
        if (self.i, self.j, self.k) == (0.0, 0.0, 0.0):
            return True

        return False

    def is_vector(self) -> bool:
        """
        Return ``True`` if the scalar part is zero and at least one of
        the vector components is nonzero. Otherwise, return ``False``.
        """
        if self.real == 0.0 and (
                self.i != 0.0 or self.j != 0.0 or self.k != 0.0):
            return True

        return False

    ## Attributes ##
    ################
    @property
    def angle(self) -> float:
        """The angle of the quaternion in radians."""
        return _atan2(self.vector.__abs__(), self.real)

    @property
    def angle_in_degrees(self) -> float:
        """The angle of the quaternion in degrees."""
        return self.angle * 180 / _pi

    angle_in_radians = angle

    @property
    def components(self) -> Tuple[float]:
        """
        The components of the quaternion as a tuple in the order
        ``(real, i, j, k)``.
        """
        return (self.real, self.i, self.j, self.k)

    @property
    def norm(self) -> float:
        """
        The norm (magnitude) of the quaternion.
        """
        return self.__abs__()

    @property
    def scalar(self) -> Quaternion:
        """The real part of the quaternion."""
        return Quaternion(self.real, 0, 0, 0)

    @property
    def vector(self) -> Quaternion:
        """The vector part of the quaternion."""
        return Quaternion(0, self.i, self.j, self.k)

    @property
    def vector_norm(self) -> float:
        """The norm of the vector part of the quaternion."""
        return self.vector.__abs__()

    @property
    def versor(self) -> Quaternion:
        """The quaternion normalized to a magnitude of one (1)."""
        versor = self
        for _ in range(10):  # Prevents an infinite loop.
            # Doesn't always divide to norm 1 on the first division.
            versor = versor/versor.__abs__()
            if versor.norm == 1.0:
                break
        return versor

    ## Class methods ##  # Only one so far.
    ###################
    @classmethod
    def from_angle(
            cls, angle: float, vector: Iterable[float],
            norm: float = 1.0, degrees: bool = True) -> Quaternion:
        """
        Return a quaternion from an angle and vector.

        Quaternions can be expressed as ``norm*(cos(theta) + u*sin(theta))``,
        where ``u`` is a 3D unit vector. This function takes an angle and a vector to
        create a quaternion. If you want a quaternion with a different norm than
        one (1), you can change the ``norm`` argument. By default, angles are entered
        in degrees. If you want to enter an angle in radians, set ``degrees`` to
        False.
        """
        if degrees:
            angle %= 360
            angle = angle * _pi / 180
        else:
            angle %= (_pi*2)

        if angle == 0.0:
            return cls(1.0)
        else:
            i, j, k = _makeListLen3(vector)
            u = cls(0, i, j, k)
            if abs(u) == 0:
                raise ValueError(
                    'The vector must be a nonzero vector.')
            elif abs(u) != 1.0:
                u = u.versor

            q = _cos(angle) + _sin(angle)*u

            if norm != 1.0:
                q = norm * q.versor

            return q
