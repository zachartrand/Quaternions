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

import sys
from numbers import Number as _Number
import math as _math
from typing import Tuple, Iterable, Iterator

from quaternions._misc import (
    deg_cos as _deg_cos,
    deg_sin as _deg_sin,
    makeListLen3 as _makeListLen3,
)

# Constants #
_LN2 = 0.69314718055994530942  # Natural log of 2.
_LARGE_FLOAT = 0.25*sys.float_info.max
_SQRT_LARGE_FLOAT = _math.sqrt(_LARGE_FLOAT)
_SMALL_FLOAT = sys.float_info.min
_MANT_DIG = sys.float_info.mant_dig
_SCALE_UP_MAX = _math.ldexp(sys.float_info.max, -_MANT_DIG)


class Quaternion():
    """
    Quaternions are an expansion of the complex numbers, where there
    are four (4) components--the real component, also known as the
    scalar part, and the imaginary components, which together are
    known as the vector part. The vector part is made up of three (3)
    components whose unit values are `i`, `j`, and `k`. The rules for
    these values are as follows:

      :math:`i^2 = j^2 = k^2 = -1`

      :math:`jk = -kj = i`

      :math:`ki = -ik = j`

      :math:`ij = -ji = k`,

    which leads to the following statement:

      :math:`ijk = -1`.

    The descriptions will reference a quaternion of the form
    :math:`a + bi + cj + dk`, where :math:`a`, :math:`b`, :math:`c`,
    and :math:`d` are real numbers.

    Parameters:
        real_component: The real component (:math:`a`) of the quaternion.
        i_component: The i component (:math:`b`) of the quaternion.
        j_component: The j component (:math:`c`) of the quaternion.
        k_component: The k component (:math:`d`) of the quaternion.

    Each component can be returned by calling the attribute of the
    same name.

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

    __slots__ = ("_real", "_i", "_j", "_k")

    ## Initialization ##
    ####################
    # Quaternions are a type of number, so make it immutable like
    # other numeric types. Use __new__ instead of __init__.
    def __new__(cls, real_component: float = 0.0, i_component: float = 0.0,
                j_component: float = 0.0, k_component: float = 0.0):
        self = object.__new__(cls)
        self._real: float = float(real_component)
        self._i: float = float(i_component)
        self._j: float = float(j_component)
        self._k: float = float(k_component)

        return self

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

        Return a complex number of the form 
        
            (self.real + self.vector_norm*1j)
        """
        imag = self.vector_norm
        
        return complex(self.real, imag)

    def __float__(self) -> float:
        """
        Return float(self).

        If the quaternion is scalar, return the scalar component.
        Otherwise, raise a ValueError.
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
        # The string returned by this function is a string that, when
        # copied to the Python interpreter, will create an identical
        # Quaternion object.
        #
        # Example: for the quaternion (1 - 2i - 3j + 4k), assuming the
        # Quaternion object is imported as Quaternion, this function
        # will return the following string: "Quaternion(1.0, -2.0, -3.0, 4.0)".

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
            """
            Return '+' or '-' based on whether x is positive or
            negative.
            """
            if _math.copysign(1.0, x) == -1.0:
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

        j_str = "".join([_sign(self.j), " ", _num_to_str(self.j.__abs__()), "j"])
        k_str = "".join([_sign(self.k), " ", _num_to_str(self.k.__abs__()), "k"])
        if self.real:
            i_str = "".join([_sign(self.i), " ", _num_to_str(self.i.__abs__()), "i"])
            q_str = " ".join([_num_to_str(self.real), i_str, j_str, k_str])
        else:
            i_str = "".join([_num_to_str(self.i), "i"])
            q_str = " ".join([i_str, j_str, k_str])

        return f"({q_str})"

    ## Boolean magic methods ##
    ###########################
    def __bool__(self) -> bool:
        """self != 0"""
        return (self.real, self.i, self.j, self.k) != (0.0, 0.0, 0.0, 0.0)

    def __eq__(self, other: Quaternion or float or complex) -> bool:
        """Return self == other."""
        # When comparing to an integer or floating-point number, the
        # real component of the quaternion must equal the int or float
        # with the vector components all equal to 0.
        if isinstance(other, (int, float)):
            return (self.real, self.i, self.j, self.k) == (other, 0.0, 0.0, 0.0)
        
        # A quaternion can only be equal to a complex number iff no
        # more than one of the vector components is nonzero. Then, the
        # real components must be equal and the nonzero vector
        # component must be equal to the imaginary component.
        elif isinstance(other, complex):
            return (self.real, self.get_imag()) == (
                other.real, other.imag)
        
        # Quaternions are equal to each other if each component is
        # equal to the other's corresponding component.
        elif isinstance(other, self.__class__):
            return (self.real, self.i, self.j, self.k) == (
                other.real, other.i, other.j, other.k)
        
        return False

    # Because several quaternions can be the same magnitude, and it's
    # too arbitrary to pick how to order them, these Boolean functions
    # are NotImplemented.
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
        # Catch non-finite components.
        if any([not _math.isfinite(component) for component in self.__iter__()]):
            # If at least one component is a NaN, return NaN.
            if any([_math.isnan(component) for component in self.__iter__()]):
                return _math.nan
            
            # If any component is an infinity and there are no NaNs, 
            # the norm is infinity.
            return _math.inf

        norm = _math.hypot(self.real, self.i, self.j, self.k)
        # Catch overflow from finite components.
        if not _math.isfinite(norm):
            raise OverflowError("absolute value too large")

        return norm

    def __ceil__(self) -> Quaternion:
        """
        Return math.ceil(self).

        Return the quaternion with all of its components rounded up
        to the nearest integer.
        """
        return self.__class__(
            _math.ceil(self.real), _math.ceil(self.i), _math.ceil(self.j), _math.ceil(self.k))

    def __floor__(self) -> Quaternion:
        """
        Return math.floor(self).

        Return the quaternion with all of its components rounded down
        to the nearest integer.
        """
        return self.__class__(
            _math.floor(self.real), _math.floor(self.i), _math.floor(self.j), _math.floor(self.k))

    def __pos__(self) -> Quaternion:
        """Return +self."""
        return self

    def __neg__(self) -> Quaternion:
        """Return -self."""
        return self.__class__(-self.real, -self.i, -self.j, -self.k)

    def __round__(self, ndigits: None or int = None) -> Quaternion:
        """
        Return round(self, ndigits).

        Rounds each component of the quaternion to the nearest integer.
        Each component is returned as a float.
        """
        return self.__class__(
            round(self.real, ndigits), round(self.i, ndigits),
            round(self.j, ndigits), round(self.k, ndigits))

    ## Mathematical magic methods with two mandatory inputs. ##
    ###########################################################
    def __add__(self, other: Quaternion or float) -> Quaternion:
        """
        Return self + other.

        If other is an int or float, it gets added to the scalar part
        of the quaternion. If other is a quaternion, the quaternions
        are added component-by-component.
        """
        # When adding to an integer or floating-point number, the
        # int/float gets added to the real component of the quaternion,
        # with the vector components remaining constant.
        if isinstance(other, (int, float)):
            if other == 0.0:
                return self
            
            real = self.real + other
            i, j, k = self.i, self.j, self.k
        # Because it is unclear which vector component the imaginary
        # part of a complex number represents, raise a TypeError
        # telling the user to convert the complex number to a
        # quaternion.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot add a {self_type} to a complex number. Make "
                + f"the complex number a {self_type} before adding.")
        # Quaternions add to each other component-by-component.
        elif isinstance(other, self.__class__):
            real = self.real + other.real
            i = self.i + other.i
            j = self.j + other.j
            k = self.k + other.k
        else:
            return NotImplemented

        return self.__class__(real, i, j, k)

    def __radd__(self, other: Quaternion or float) -> Quaternion:
        """Return other + self."""
        # Addition with quaternions is commutative, so we can use the
        # normal __add__ method for ints and floats.
        if isinstance(other, (int, float)):
            return self.__add__(other)
        
        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot add a complex number to a {self_type}. Make "
                + f"the complex number a {self_type} before adding.")
        
        return NotImplemented

    def __sub__(self, other: Quaternion or float) -> Quaternion:
        """
        Return self - other.

        If other is an int or float, it is subtracted from the scalar
        part of the quaternion. If other is a quaternion, the second
        quaternion is subtracted by the first component-by-component.
        """
        # Like addition, if an int/float is subtracted from a
        # quaternion, it gets subtracted from the real component with
        # the other components remaining constant.
        if isinstance(other, (int, float)):
            if other == 0.0:
                return self
            
            real = self.real - other
            i, j, k = self.i, self.j, self.k
        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot subtract a complex number from a {self_type}. Make "
                + f"the complex number a {self_type} before subtracting.")
        
        # A quaternion subtracted from another quaternion subtracts
        # component by component.
        elif isinstance(other, self.__class__):
            real = self.real - other.real
            i = self.i - other.i
            j = self.j - other.j
            k = self.k - other.k
        else:
            return NotImplemented
        
        return self.__class__(real, i, j, k)

    def __rsub__(self, other: Quaternion or float) -> Quaternion:
        """
        Return other - self.

        If other is an int or float, the scalar part of the quaternion
        is subtracted from other, and the signs of the vector
        components are switched.
        """
        # When a quaterion is subtracted from an int/float, the real
        # component is subtracted from the int/float, and the other
        # components are negated.
        if isinstance(other, (int, float)):
            return -self.__add__(other)

        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot subtract a {self_type} from a complex number. Make "
                + f"the complex number a {self_type} before subtracting.")

        return NotImplemented

    def __mul__(self, other: Quaternion or float) -> Quaternion:
        """
        Return self * other.

        If the Quaternion is multiplied by an int or float, the
        int/float is distributed to each component. If other is a
        quaternion, this multiplies following the rules of quaternion
        multiplication. Multiplication is not commutative; q1 * q2
        usually does not have the same result as q2 * q1.
        """
        # When a quaternion is multiplied by an int/float, that number
        # is distributed to each component.
        if isinstance(other, (int, float)):
            if other == 1.0:
                return self
            
            real, i, j, k = (
                other*self.real, other*self.i, other*self.j, other*self.k)
        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot multiply a {self_type} by a complex number. Make the "
                + f"complex number a {self_type} before multiplying.")
        #
        # For quaternion multiplication with other quaternions, order
        # matters. The resulting sum of products is based on how the
        # unit vector components multiply with each other. See the
        # docstring of the Quaternion class to see how these values
        # multiply with each other.
        elif isinstance(other, self.__class__):
            max_component = max(*self.abs_components(), *other.abs_components())
            if max_component > _SQRT_LARGE_FLOAT:
                # Prevent potential overflow by scaling down values 
                # by largest component, calculating the sums, and then 
                # rescaling back up.
                real = max_component*(
                    (max(self.real, other.real)/max_component)*min(self.real, other.real)
                    - (max(self.i, other.i)/max_component)*min(self.i, other.i)
                    - (max(self.j, other.j)/max_component)*min(self.j, other.j)
                    - (max(self.k, other.k)/max_component)*min(self.k, other.k)
                )
                i = max_component*(
                    (max(self.real, other.i)/max_component)*min(self.real, other.i)
                    + (max(self.i, other.real)/max_component)*min(self.i, other.real) 
                    + (max(self.j, other.k)/max_component)*min(self.j, other.k)
                    - (max(self.k, other.j)/max_component)*min(self.k, other.j)
                )
                j = max_component*(
                    (max(self.real, other.j)/max_component)*min(self.real, other.j)
                    + (max(self.j, other.real)/max_component)*min(self.j, other.real) 
                    + (max(self.k, other.i)/max_component)*min(self.k, other.i)
                    - (max(self.i, other.k)/max_component)*min(self.i, other.k)
                )
                k = max_component*(
                    (max(self.real, other.k)/max_component)*min(self.real, other.k)
                    + (max(self.k, other.real)/max_component)*min(self.k, other.real)
                    + (max(self.i, other.j)/max_component)*min(self.i, other.j)
                    - (max(self.j, other.i)/max_component)*min(self.j, other.i)
                )
            
            else:
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
        
        return self.__class__(real, i, j, k)

    def __rmul__(self, other: Quaternion or float) -> Quaternion:
        """Return other * self."""
        if isinstance(other, (int, float)):
            return self.__mul__(other)  # Scalar multiplication commutes.
        
        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot multiply a complex number by a {self_type}. Make "
                + f"the complex number a {self_type} before multiplying.")
        
        return NotImplemented

    def __floordiv__(self, other: float) -> Quaternion:
        """
        Return self // other.

        If other is an int or float, return the quaternion with each
        component floor divided by other.
        """
        if isinstance(other, (int, float)):
            if other == 1.0:
                return self
            
            return self.__class__(
                self.real // other, self.i // other, self.j // other,
                self.k // other)

        return NotImplemented

    def __truediv__(self, other: Quaternion or float) -> Quaternion:
        """
        Return self/other.

        For division q1/q2, this assumes the order of multiplication
        is q1 * 1/q2. To left-multiply the denominator, enter 1/q2 * q1
        or q2.inverse() * q1.
        """
        # Division by an int/float distibutes through each component
        # just like scalar multiplication.
        if isinstance(other, (int, float)):
            if other == 1.0:
                return self
            
            inv_other = other**-1
            return self.__class__(
                self.real*inv_other, self.i*inv_other,
                self.j*inv_other, self.k*inv_other)

        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot divide a {self_type} by a complex number. Make "
                + f"the complex number a {self_type} before dividing.")

        # Multiply by the multiplicative inverse of the second
        # quaternion. See the inverse() method to see how this is
        # calculated.
        elif isinstance(other, self.__class__):
            if other.is_scalar():
                if other.real == 0.0:  # other is the zero quaternion.
                    raise ZeroDivisionError("Quaternion division by zero.")
                
                return self.__truediv__(other.real)
            
            elif self.is_scalar():
                return other.inverse().__mul__(self.real)
            
            return self.__mul__(other.inverse())

        return NotImplemented

    def __rtruediv__(self, other: Quaternion or float) -> Quaternion:
        """Return other/self."""
        # Multiply the multiplicative inverse of the quaternion by the
        # int/float.
        if isinstance(other, (int, float)):
            return self.inverse().__mul__(other)

        # See __add__ for complex logic.
        elif isinstance(other, complex):
            self_type = self.__class__.__qualname__
            raise TypeError(
                f"Cannot divide a complex number by a {self_type}. Make "
                + f"the complex number a {self_type} before dividing.")
        
        return NotImplemented
    
    def __pow__(self, other: Quaternion or float,
                mod: None = None) -> Quaternion:
        """
        Return pow(self, other, mod).

        Return self to the power of other, where other is a real
        exponent.
        """
        # Modular arithmetic doesn't make sense with quaternions.
        # Raises a ValueError.
        if mod is not None:
            raise ValueError(
                f"{self.__class__.__qualname__}s are not set up "
                + "for modular arithmetic.")
        #
        # The code where `other` is an int or float has been updated. 
        # The old code ran on O(other) time, while the new code runs 
        # on O(log(other)) time. See version 1.1.3 for the old code.
        one = self.__class__(1, 0, 0, 0)
        if (isinstance(other, int)
                or (isinstance(other, float) and other.is_integer())):
            other = int(other)
            if other == 0:
                return one
            
            # Avoid quaternion math if the quaternion is a scalar.
            if self.is_scalar():
                # Special cases for the zero quaternion.
                if self.real == 0.0:
                    if other > 0:
                        return self
                    elif other < 0:
                        raise ZeroDivisionError(
                            f"0.0 cannot be raised to a negative power.")
                
                return self.__class__(pow(self.real, other), 0, 0, 0)
            
            if other < 0:
                base = self.inverse()
                other = abs(other)
            else:
                base = self
            
            result = one
            while other > 1:
                if other % 2 == 1:
                    result = result*base
                
                base = base.squared()
                other = other//2
            
            result = result*base
            
            return result
        
        elif isinstance(other, float):
            # Avoid quaternion math if the quaternion is a scalar.
            if self.is_scalar():
                # Special cases for the zero quaternion.
                if self.real == 0.0:
                    if other > 0.0:
                        return self
                    elif other < 0.0:
                        raise ZeroDivisionError(
                            f"0.0 cannot be raised to a negative power.")
                
                return self.__class__(pow(self.real, other), 0, 0, 0)
            
            other_int = int(other)
            other = other - other_int

            result = self.__pow__(other_int)
            theta = other * self.angle
            new_norm = pow(self.__abs__(), other)
            v = self.vector
            v_norm = self.vector.norm
            result = result*(_math.cos(theta)*new_norm + v*_math.sin(theta)*(new_norm/v_norm))
            
            return result
        
        elif isinstance(other, self.__class__):
            # Avoid more complicated math if the exponent is Real.
            if other.is_scalar():
                return self.__pow__(other.real)
            
            # A quaternion (q) raised to a quaternion power (p) is
            # given by the following formula:
            #
            #     q**p = exp(ln(q)*p),
            #
            # where ln is the natural logarithm, ln(q) = ln(||q||) +
            # v*phi, and exp is the exponential function (see the exp()
            # function in the qmath module).
            elif self.is_scalar():
                # Copy the functionality of complex numbers 
                # for zero raised to a quaternion:
                if self.real == 0.0:  # self is the zero quaternion.
                    raise ZeroDivisionError("Zero cannot be raised to a quaternion.")
                
                q = _math.log(abs(self.real))*other
            else:
                v_norm = self.vector_norm
                ln = self.log_norm()*v_norm + self.vector*self.angle
                q = ln * (other/v_norm)
            
            q_vector = q.vector
            q_vnorm = q.vector_norm
            pow_q = (
                (_math.exp(q.real)/q_vnorm) * (q_vector.cos_norm()*q_vnorm + q_vector*q_vector.sin_norm()))
            
            return pow_q
        
        return NotImplemented
    
    def __rpow__(self, other, mod=None):
        """Return other**self"""
        if isinstance(other, (int, float)):
            return Quaternion(other, 0, 0, 0).__pow__(self, mod)
        
        return NotImplemented

    # Modular arithmetic does not make sense for quaternions, so the
    # __mod__() method is NotImplemented.
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
        if self.is_scalar():
            return self
        
        return self.__class__(self.real, -self.i, -self.j, -self.k)

    def get_imag(self) -> float or None:
        """
        Return the imaginary component of the quaternion if only one
        of the imaginary components is nonzero. If the quaternion is
        scalar, return ``0.0``. Otherwise, return ``None``.
        """
        if self.is_scalar():
            return 0.0
        
        elif self.is_complex():
            for component in (self.i, self.j, self.k):
                if component != 0.0:
                    return component
        
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
        # ==============================================================
        # Added an algorithm similar to Python's algorithm for complex
        # numbers to avoid overflow and underflow errors. Using this
        # algorithm avoids squaring numbers and scales down each 
        # component by the largest of the components, making the base 
        # numbers to work with between 0 and 1 in magnitude.
        #
        # UPDATE 20 Jul 2022: Changed float divisions to multiplications
        # for faster calculations
        # ==============================================================
        #
        # Check for special cases of norm zero (0.0) and norm one (1.0).
        if self.is_scalar():
            if self.real == 0.0:  # self is the zero quaternion.
                raise ZeroDivisionError("The zero quaternion has no defined inverse.")
            
            return self.__class__(pow(self.real, -1), 0, 0, 0)
        
        elif self.norm == 1.0:
            return self.conjugate()
        
        inv_max_component = (max(self.abs_components()))**-1
        real_ratio = self.real * inv_max_component
        i_ratio = self.i * inv_max_component
        j_ratio = self.j * inv_max_component
        k_ratio = self.k * inv_max_component
        inv_denom = (real_ratio*self.real + i_ratio*self.i
                    + j_ratio*self.j + k_ratio*self.k)**-1

        q_inverse = (
            self.__class__(real_ratio*inv_denom, -i_ratio*inv_denom,
                        -j_ratio*inv_denom, -k_ratio*inv_denom))

        return q_inverse
    
    def squared(self):
        """Return self**2."""
        max_component = max(self.abs_components())
        inv_max_component = 1/max_component
        real = _math.fsum([
            (self.real*inv_max_component)*self.real, 
            (-self.i*inv_max_component)*self.i, 
            (-self.j*inv_max_component)*self.j, 
            (-self.k*inv_max_component)*self.k
        ])
        real = real*max_component

        i = 2.0*self.real*self.i
        j = 2.0*self.real*self.j
        k = 2.0*self.real*self.k

        return self.__class__(real, i, j, k)

    def unit_quaternion(self) -> Quaternion:
        """
        Return the quaternion normalized to magnitude one (1).

        If the quaternion is a zero (0) quaternion, return the zero
        quaternion.
        """
        if self.is_zero():
            return self

        return self.versor

    def unit_vector(self) -> Quaternion:
        """
        Return the vector part of the quaternion normalized to a
        magnitude of one (1.0). Return the zero quaternion if the
        magnitude of the quaternion is zero (0.0).
        """
        vector = self.vector
        if vector.is_zero():
            return self
        
        return vector.versor
    
    def log_norm(self) -> float:
        """
        Return the natural logarithm of the norm of self.

        This tends to be more accurate than 
        
        >>> math.log(self.norm)
        """
        # This function takes the size of the absolute values of the 
        # components into consideration when calculating the logarithm
        # of the norm. If any of the components are large enough to 
        # cause a potential overflow, small enough to cause potential 
        # underflow, or the norm itself is close to one (1.0), each 
        # case has its own algorithm to account for errors.
        #
        # Because log1p() is more accurate than log() for inputs close to
        # 1.0, when the norm of the quaternion is between 0.71 and 1.73
        # inclusive, log1p is automatically used. This range is used based
        # on the cmath module, which uses log1p for complex numbers with
        # magnitudes in this range.
        abs_components = self.abs_components()
        max_component = max(abs_components)

        # Check if any of the components could cause an overflow.
        if max_component > _LARGE_FLOAT:
            result = _math.log(_math.hypot(*[0.5*x for x in abs_components])) + _LN2
        elif max_component < _SMALL_FLOAT:
            # Check if any components are subnormal.
            if max_component > 0.0:
                result = (
                    _math.log(_math.hypot(*[_math.ldexp(x, _MANT_DIG) for x in abs_components]))
                    - _MANT_DIG*_LN2
                )
            else:  # norm is zero (0.0).
                return _math.log(0.0)
        else:
            norm = self.__abs__()
            if 0.71 <= norm and norm <= 1.73:
                # This takes advantage of the difference of perfect
                # squares to avoid potential loss of precision in
                # calculating the norm - 1 of the quaternion. This
                # formula is taken from the cmath log function,
                # modified for four (4) components.
                abs_components.remove(max_component)
                result = 0.5*_math.log1p(
                    (max_component - 1)*(max_component + 1)
                    + abs_components[0]*abs_components[0]
                    + abs_components[1]*abs_components[1]
                    + abs_components[2]*abs_components[2]
                )
            else:
                result = _math.log(norm)

        return result
    
    def from_complex(self, z: complex) -> Quaternion:
        """
        Return a Quaternion from a complex number and the vector of self.
        
        If ``u == self.unit_vector()``, this is equivalent to 
        
            ``Quaternion(z.real, z.imag*u.i, z.imag*u.j, z.imag*u.k)``
        """
        vector = self.get_vector_components()
        v_norm = self.vector_norm

        return Quaternion(z.real, *[x if x == 0.0 else x*(z.imag/v_norm) for x in vector])
    
    def abs_components(self):
        """
        Return a list of the absolute values of the components of self.
        """
        return [abs(x) for x in self.__iter__()]
    
    def cos_norm(self) -> float:
        """Return the cosine of the norm of self."""
        max_component = max(self.abs_components())
        if max_component > _LARGE_FLOAT:
            # Calculate half the vector norm and use 
            # the double angle formula for cosine.
            half_norm = _math.hypot(0.5*self.real, 0.5*self.i, 0.5*self.j, 0.5*self.k)
            cosine = _math.cos(half_norm)
            sine = _math.sin(half_norm)

            return cosine*cosine - sine*sine
        
        return _math.cos(self.__abs__())
    
    def sin_norm(self) -> float:
        """Return the sine of the norm of self."""
        max_component = max(self.abs_components())
        if max_component > _LARGE_FLOAT:
            # Calculate half the vector norm and use 
            # the double angle formula for sine.
            half_norm = _math.hypot(0.5*self.real, 0.5*self.i, 0.5*self.j, 0.5*self.k)

            return 2.0*_math.sin(half_norm)*_math.cos(half_norm)
        
        return _math.sin(self.__abs__())

    ## Boolean methods ##
    #####################
    def is_complex(self) -> bool:
        """
        Return ``True`` if only one of the *i*, *j*, and *k* components
        is nonzero. Otherwise, return ``False``.
        """
        if self.is_scalar():
            return False
        
        return (0.0, 0.0) in (
                    (self.i, self.j), (self.j, self.k), (self.i, self.k))

    def is_scalar(self) -> bool:
        """
        Return ``True`` if the vector components all equal zero.
        Otherwise, return ``False``.
        """
        return (self.i, self.j, self.k) == (0.0, 0.0, 0.0)

    def is_vector(self) -> bool:
        """
        Return ``True`` if the scalar part is zero and at least one of
        the vector components is nonzero. Otherwise, return ``False``.
        """
        return ((self.real == 0.0) and (
                (self.i != 0.0) or (self.j != 0.0) or (self.k != 0.0)))
    
    def is_zero(self) -> bool:
        return not self.__bool__()
    
    def is_not_zero(self) -> bool:
        return self.__bool__()

    ## Attributes ##
    ################
    @property
    def angle(self) -> float:
        """The angle of the quaternion in radians."""
        # NOTE: This uses atan2 over acos due to inaccurate calculations
        # near 0 and pi radians for acos. See here:
        #
        #   https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Numerical_accuracy
        #
        # This makes the formula for the angle
        # 
        #     atan2(self.vector_norm, self.real)
        #
        # The main issue with this formula is that the vector norm, 
        # being a calculation of the hypotenuse of the vector 
        # components, can overflow/underflow with large/small absolute 
        # values in the vector components. Since arctangent is a 
        # ratio, large and small magnitude components are rescaled to 
        # get accurate values for the quaternion's angle.
        real = self.real
        vector_components = self.abs_components()[1:]
        max_component = max(vector_components)
        if max_component > _LARGE_FLOAT:
            # Risk of overflow, reduce all components by half.
            vector_components = [0.5*x for x in vector_components]
            real = 0.5*real
        elif max_component < _SMALL_FLOAT:
            # Find subnormal floats.
            if max_component > 0.0:
                # At least one component is subnormal, so scale the 
                # components up if `real` is sufficiently low.
                if abs(real) <= _SCALE_UP_MAX:
                    vector_components = [_math.ldexp(x, _MANT_DIG) for x in vector_components]
                    real = _math.ldexp(real, _MANT_DIG)
                
            else:  # self is scalar.
                return _math.atan2(0.0, real)
        
        vector_norm = _math.hypot(*vector_components)

        return _math.atan2(vector_norm, real)

    @property
    def angle_in_degrees(self) -> float:
        """The angle of the quaternion in degrees."""
        return _math.degrees(self.angle)

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
        """The norm (magnitude) of the quaternion."""
        return self.__abs__()

    @property
    def scalar(self) -> Quaternion:
        """The real part of the quaternion."""
        return self.__class__(self.real, 0, 0, 0)

    @property
    def vector(self) -> Quaternion:
        """The vector part of the quaternion."""
        return self.__class__(0, self.i, self.j, self.k)

    @property
    def vector_norm(self) -> float:
        """The norm of the vector part of the quaternion."""
        # Be sure to go through __abs__ method to ensure overflow is caught.
        return self.vector.__abs__()

    @property
    def versor(self) -> Quaternion:
        """The quaternion normalized to a magnitude of one (1)."""
        # Deal with NaNs.
        if any([_math.isnan(component) for component in self.__iter__()]):
            return Quaternion(*[x if x == 0.0 else _math.nan for x in self.__iter__()])
        
        # Deal with large and small norms.
        max_component = max(self.abs_components())
        if max_component > _LARGE_FLOAT:
            if _math.isinf(max_component):
                components = [_math.copysign(1.0, x) if _math.isinf(x) else _math.copysign(0.0, x) for x in self.__iter__()]
                return self.__class__(*components).versor
            else:
                versor = 0.5*self
        elif max_component < _SMALL_FLOAT:
            if max_component > 0.0:
                versor = self*pow(2, _MANT_DIG)
            else:  # self is the zero quaternion.
                raise ZeroDivisionError(
                    "The zero quaternion cannot be converted to norm one (1).")
        else:
            versor = self
        
        versor = versor / versor.__abs__()
        if versor.__abs__() != 1.0:
            # Doesn't always divide to norm 1 on the first division.
            # Attempt at most ten (10) times to reduce to norm one (1).
            scale = 1
            primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
            for prime in primes:
                # Because of the nature of floating point numbers, the 
                # first division doesn't always reduce to norm one. If 
                # the norm of versor is not equal to one after 
                # division, scale it up by a factor made up of the 
                # first odd primes.
                scale *= prime
                versor = versor*scale
                versor = versor / versor.__abs__()
                if versor.__abs__() == 1.0:
                    break
        
        return versor

    ## Class methods ##  # Only one so far.
    ###################
    @classmethod
    def from_angle(cls, angle: float, vector: Iterable[float],
                    norm: float = None, degrees: bool = True) -> Quaternion:
        """
        Return a quaternion from an angle and vector.

        Quaternions can be expressed as ``norm*(cos(theta) +
        u*sin(theta))``, where ``u`` is a 3D unit vector. This function
        takes an angle and a vector to create a quaternion. If you want
        a quaternion with a specific magnitude, you can change
        the ``norm`` argument. If no argument is given for `norm`, the 
        resulting quaternion will have a norm equal to the magnitude of 
        `vector`. By default, angles are entered in degrees. If you 
        want to enter an angle in radians, set ``degrees`` to False.
        """
        if degrees:
            cos = _deg_cos
            sin = _deg_sin
        else:
            cos = _math.cos
            sin = _math.sin
            
        i, j, k = _makeListLen3(vector)
        # If the vector is the zero vector, only the scalar part
        # remains.
        if i == 0.0 and j == 0.0 and k == 0.0:
            if norm is not None:
                return cls(norm)

            return cls(cos(angle))

        vector_norm = _math.hypot(i, j, k)
        
        sine = sin(angle)
        q = cls(cos(angle)*vector_norm, i*sine, j*sine, k*sine)

        if norm is not None:
            q = (norm/vector_norm)*q

        return q


_Number.register(Quaternion)
