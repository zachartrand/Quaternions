"""
Created on Sat May 29 23:23:00 2021

@author: Zach Chartrand <zachartrand999@gmail.com>
"""

from math import acosh as _acosh, hypot as _hypot


class HyperbolicQuaternion():
    """"""
    def __init__(self, real_comp: int or float=0, i_comp: int or float=0,
                 j_comp: int or float=0, k_comp: int or float=0):
        self.real = float(real_comp)
        self.i = float(i_comp)
        self.j = float(j_comp)
        self.k = float(k_comp)

    def __iter__(self):
        """Returns an iterator of the components of the quaternion."""
        return iter([self.real, self.i, self.j, self.k])

    def __reversed__(self):
        """
        Returns an iterator of the components of the quaternion
        in reverse order.
        """
        return reversed([self.real, self.i, self.j, self.k])

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

    def __pos__(self):
        """Return +self."""
        return self

    def __neg__(self):
        """Return -self."""
        return HyperbolicQuaternion(-self.real, -self.i, -self.j, -self.k)

    def __abs__(self):
        max_component = max(
            abs(self.real), abs(self.i), abs(self.j), abs(self.k))
        real_ratio = self.real / max_component
        i_ratio = self.i / max_component
        j_ratio = self.j / max_component
        k_ratio = self.k / max_component
        hypot = (real_ratio * self.real - i_ratio * self.i - j_ratio * self.j
            - k_ratio * self.k)

        norm = max_component**0.5 * hypot**0.5
        return norm

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

    def __bool__(self):
        """self != 0"""
        if (self.real, self.i, self.j, self.k) == (0.0, 0.0, 0.0, 0.0):
            return False

        return True

    def __eq__(self, other):
        """Return self == other."""
        if isinstance(other, (int, float)):
            return (self.real, self.i, self.j, self.k) == (other, 0.0, 0.0, 0.0)
        elif isinstance(other, HyperbolicQuaternion):
            return (self.real, self.i, self.j, self.k) == (
                other.real, other.i, other.j, other.k)

        return False

    def __hash__(self):
        """Return hash(self)."""
        if self.is_scalar():
            return hash(self.real)
        else:
            return hash((self.real, self.i, self.j, self.k, 'hyperbolic'))

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
        elif isinstance(other, Quaternion):
            real = self.real + other.real
            i = self.i + other.i
            j = self.j + other.j
            k = self.k + other.k
        else:
            return NotImplemented

        return HyperbolicQuaternion(real, i, j, k)

    def __radd__(self, other):
        """Return other+self."""
        if isinstance(other, (int, float)):
            return self.__add__(other)

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
        elif isinstance(other, Quaternion):
            real = self.real - other.real
            i = self.i - other.i
            j = self.j - other.j
            k = self.k - other.k
        else:
            return NotImplemented

        return HyperbolicQuaternion(real, i, j, k)

    def __rsub__(self, other):
        """
        Return other-self.

        If other is an int or float, the scalar part of the quaternion
        is subtracted from other, and the signs of the vector
        components are switched.
        """
        if isinstance(other, (int, float)):
            real = other - self.real
            return HyperbolicQuaternion(real, -self.i, -self.j, -self.k)

        return NotImplemented

    def __mul__(self, other):
        """Return self * other."""
        if isinstance(other, HyperbolicQuaternion):
            real = (self.real*other.real + self.i*other.i + self.j*other.j
            + self.k*other.k)
            i = (self.real*other.i + self.i*other.real + self.j*other.k
                 - self.k*other.j)
            j = (self.real*other.j + self.j*other.i + self.k*other.i
                 - self.i*other.k)
            k = (self.real*other.k + self.k*other.real + self.i*other.j
                 - self.j*other.i)

            return HyperbolicQuaternion(real, i, j, k)

        elif isinstance(other, (int, float)):
            return HyperbolicQuaternion(self.real*other, self.i*other,
                self.j* other, self.k*other)

        return NotImplemented

    def __rmul__(self, other):
        """Return other * self."""
        if isinstance(other, (int, float)):
            return self.__mul__(other)

    def __truediv__(self, other):
        """Return self/other."""
        if isinstance(other, (int, float)):
            return HyperbolicQuaternion(self.real/other, self.i/other,
                self.j/other, self.k/other)
        elif isinstance(other, HyperbolicQuaternion):
            return self.__mul__(other.inverse())

    def __rtruediv__(self, other):
        """Return other/self."""
        if isinstance(other, (int, float)):
            return self.inverse().__mul__(other)

        return NotImplemented

    @property
    def scalar(self) -> float:
        """Return the real part of the quaternion."""
        return self.real

    @property
    def vector(self):
        """Return the vector part of the quaternion."""
        return HyperbolicQuaternion(0, self.i, self.j, self.k)

    @property
    def norm(self):
        """Returns the norm (magnitude) of the quaternion."""
        return abs(self)

    @property
    def v_magnitude(self):
        """Returns the magnitude of the vector part of the quaternion."""
        return _hypot(self.i, self.j, self.k)

    @property
    def angle(self):
        """
        Returns the hyperbolic angle of the quaternion.
        """
        return _acosh(self.versor.real)

    @property
    def versor(self):
        """Returns the quaternion normalized to a magnitude of one (1)."""
        versor = self
        while abs(versor) != 1.0:  # Doesn't always divide to norm 1 on the
                                   # first division.
            versor = versor/abs(versor)
        return versor

    def inverse(self):
        """Return 1/self."""
        max_component = max(
            abs(self.real), abs(self.i), abs(self.j), abs(self.k))
        real_ratio = self.real / max_component
        i_ratio = self.i / max_component
        j_ratio = self.j / max_component
        k_ratio = self.k / max_component
        denom = (real_ratio * self.real - i_ratio * self.i - j_ratio * self.j
            - k_ratio * self.k)

        hq_inverse = HyperbolicQuaternion(
            real_ratio, -i_ratio, -j_ratio, -k_ratio) / denom

        return hq_inverse

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
        if self.real == 0 and (self.i != 0.0 or self.j != 0.0 or self.k != 0.0):
            return True

        return False
