"""
Unittest module for the Quaternion class.
"""

import sys
import unittest

from math import hypot, log, pi, cos, sin
from random import random

from quaternions import Quaternion

large_float = sys.float_info.max
small_float = sys.float_info.min
inf = float("inf")

q1 = Quaternion(1, -2, -3, 4)
q2 = Quaternion(1, 4, -3, -2)

q1_plus_q2 = Quaternion(2, 2, -6, 2)
q1_minus_q2 = Quaternion(0, -6, 0, 6)
q2_minus_q1 = Quaternion(0, 6, 0, -6)
q1_times_q2 = Quaternion(8.0, 20.0, 6.0, 20.0)
q2_times_q1 = Quaternion(8.0, -16.0, -18.0, -16.0)
q1_dividedby_q2 = Quaternion(-0.2, -0.8, -0.4, -0.4)
q2_inverse_times_q1 = Quaternion(-0.2, 0.4, 0.4, 0.8)
q2_dividedby_q1 = Quaternion(-0.2, 0.8, 0.4, 0.4)
q1_inverse_times_q2 = Quaternion(-0.2, -0.4, -0.4, -0.8)

zero = Quaternion(0)
one = Quaternion(1)
largest_q = large_float*Quaternion(1, 1, 1, 1)
coefficients = [random() for _ in range(4)]
random_q = Quaternion(*coefficients)
large_q = large_float*random_q
small_q = small_float*random_q
inf_q = Quaternion(inf, large_float, -inf, 0.5*small_float)


class TestQuaternions(unittest.TestCase):

    def test_abs(self):
        abs_q1and2 = hypot(1, 2, 3, 4)

        with self.assertRaises(OverflowError):
            abs(largest_q)
        self.assertEqual(abs(one), 1.0)
        self.assertEqual(abs(zero), 0.0)
        self.assertEqual(abs(q1), abs_q1and2)
        self.assertEqual(abs(q2), abs_q1and2)
        self.assertAlmostEqual(abs(large_float*random_q.versor), large_float, delta=large_float*1e-15)
        self.assertAlmostEqual(abs(small_float*random_q.versor), small_float, delta=small_float*1e-15)

    
    def test_addition(self):
        self.assertEqual(q1+q2, q1_plus_q2)
        self.assertEqual(q2+q1, q1_plus_q2)

    def test_subtraction(self):
        self.assertEqual(q1-q2, q1_minus_q2)
        self.assertEqual(q2-q1, q2_minus_q1)

    def test_multiplication(self):
        self.assertEqual(q1*q2, q1_times_q2)
        self.assertEqual(q2*q1, q2_times_q1)
        self.assertAlmostEqual(1/q2*q1, q2_inverse_times_q1, places=12)
        self.assertAlmostEqual(1/q1*q2, q1_inverse_times_q2, places=12)
        self.assertEqual(q1*q1.inverse(), one)
        self.assertEqual(q1.inverse()*q1, one)
        self.assertEqual(q2*q2.inverse(), one)
        self.assertEqual(q2.inverse()*q2, one)
        # Multiply by zero.
        self.assertEqual(zero*q1, zero)
        self.assertEqual(q1*zero, zero)
        self.assertEqual(0.0*q1, zero)
        # Multiply by one.
        self.assertIs(1*q1, q1)
        self.assertIs(q2*1.0, q2)
        self.assertEqual(q1*one, q1)
        self.assertEqual(one*q2, q2)

    def test_division(self):
        self.assertAlmostEqual(q1/q2, q1_dividedby_q2, places=12)
        self.assertAlmostEqual(q2/q1, q2_dividedby_q1, places=12)
        self.assertEqual(q1//2, Quaternion(0, -1, -2, 2))
        # Divide by zero.
        with self.assertRaises(ZeroDivisionError):
            q1/zero
        with self.assertRaises(ZeroDivisionError):
            q2/0.0

    def test_power(self):
        # Quaternion squared.
        self.assertEqual(q1*q1, q1**2)
        # Test pow(q, -1) == q.inverse()
        self.assertEqual(q1.inverse(), q1**-1)
        # q == pow(q, 1)
        self.assertEqual(q1, q1**1)
        # q**0 == 1
        self.assertEqual(q1**0, one)
        self.assertEqual(q1**zero, one)
        # pow(q, -2) == q.inverse()**2
        self.assertAlmostEqual(q1**-2, q1.inverse()*q1.inverse(), places=17)
        # Zero quaternion test cases.
        self.assertEqual(pow(zero, zero), one)
        self.assertEqual(pow(zero, one), zero)
        with self.assertRaises(ZeroDivisionError):
            pow(zero, q1)
        with self.assertRaises(ZeroDivisionError):
            pow(zero, -1)
        with self.assertRaises(ZeroDivisionError):
            pow(zero, -one)
    
    def test_boolean_methods(self):
        # Test bool() function.
        self.assertTrue(bool(small_float*random_q))
        self.assertFalse(bool(zero))
        self.assertTrue(bool(one))
        # Test is_scalar() method.
        self.assertTrue(zero.is_scalar())
        self.assertTrue(one.is_scalar())
        self.assertFalse(q1.is_scalar())

    def test_is_zero(self):
        self.assertTrue(zero.is_zero())
        self.assertFalse(one.is_zero())
        # Test is_not_zero() method.
        self.assertFalse(zero.is_not_zero())
        self.assertTrue(one.is_not_zero())
        
    def test_is_vector(self):
        # Test is_vector() method.
        self.assertTrue(Quaternion(0, 1, 2, 3).is_vector())
        self.assertFalse(zero.is_vector())
        self.assertFalse(one.is_vector())
        self.assertFalse(q1.is_vector())
        self.assertTrue(Quaternion(0, 1, 0, 0).is_vector())
        self.assertTrue(Quaternion(0, 0, 1, 0).is_vector())
        self.assertTrue(Quaternion(0, 0, 0, 1).is_vector())
        self.assertTrue(Quaternion(0, 1, 1, 0).is_vector())
        self.assertTrue(Quaternion(0, 1, 0, 1).is_vector())
        self.assertTrue(Quaternion(0, 0, 1, 1).is_vector())
        self.assertFalse(Quaternion(1, 1, 0, 0).is_vector())
        self.assertFalse(Quaternion(1, 0, 0, 1).is_vector())
    
    def test_is_complex(self):
        # Test is_complex() method.
        self.assertFalse(zero.is_complex())
        self.assertFalse(one.is_complex())
        self.assertTrue(Quaternion(0, 1, 0, 0).is_complex())
        self.assertTrue(Quaternion(0, 0, 1, 0).is_complex())
        self.assertTrue(Quaternion(0, 0, 0, 1).is_complex())
        self.assertTrue(Quaternion(1, 1, 0, 0).is_complex())
        self.assertTrue(Quaternion(1, 0, 1, 0).is_complex())
        self.assertTrue(Quaternion(1, 0, 0, 1).is_complex())
        self.assertFalse(Quaternion(1, 0, 1, 1).is_complex())

    def test_from_angle(self):
        # Test one case against expected values to 6 decimal places.
        self.assertAlmostEqual(
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False),
            Quaternion(2.500000, 1.157275, 2.314550, 3.471825), places=6)
        # Test the degrees flag as both True and False and make sure
        # the result is (almost) the same.
        self.assertAlmostEqual(Quaternion.from_angle(60, [1, 2, 3], norm=5),
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False), places=15)
        # Test that an angle of zero produces a scalar quaternion.
        self.assertEqual(
            Quaternion.from_angle(0, [1, 2, 3], norm=5), Quaternion(5))
        # Test that an angle of zero produces a scalar quaternion
        # even with an irrational norm.
        self.assertEqual(
            Quaternion.from_angle(0, [1, 2, 3], norm=pi),
            pi/hypot(1, 2, 3)*Quaternion.from_angle(0, [1, 2, 3]))
        # Test that the norm parameter works properly, creating a
        # quaternion of that norm.
        self.assertAlmostEqual(
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False).norm, 5.0, places=15)
    
    def test_log_norm(self):
        self.assertAlmostEqual(largest_q.log_norm(), log(large_float) + log(2), places=15)
        self.assertAlmostEqual(large_q.log_norm(), 
            log(large_float) + log(hypot(*coefficients)), places=12)
        self.assertAlmostEqual(small_q.log_norm(), 
            log(small_float) + log(hypot(*coefficients)), places=12)
        self.assertAlmostEqual(
            Quaternion.from_angle(61, [0.9, 1, 1.1], norm=1.0001).log_norm(), 
            9.999500033335019e-05, places=15)

        with self.assertRaises(ValueError):
            zero.log_norm()
    
    def test_cos_norm(self):
        norms = [0.0, pi/6, pi/3, pi/2, pi, 1.5*pi, 2.0*pi, small_float, large_float]
        for norm in norms:
            kwargs = dict(delta = norm*10**-15)
            self.assertAlmostEqual((norm*random_q.versor).cos_norm(), cos(norm), **kwargs)
            self.assertAlmostEqual((-norm*random_q.versor).cos_norm(), cos(norm), **kwargs)
    
    def test_sin_norm(self):
        norms = [0.0, pi/6, pi/3, pi/2, pi, 1.5*pi, 2.0*pi, small_float, large_float]
        for norm in norms:
            kwargs = dict(delta = norm*10**-15)
            self.assertAlmostEqual((norm*random_q.versor).sin_norm(), sin(norm), **kwargs)
            self.assertAlmostEqual((-norm*random_q.versor).sin_norm(), sin(norm), **kwargs)

    def test_versor(self):
        self.assertAlmostEqual(largest_q.versor, Quaternion(1, 1, 1, 1).versor, places=15)
        self.assertAlmostEqual(large_q.versor, random_q.versor, places=15)
        self.assertAlmostEqual(small_q.versor, random_q.versor, places=15)
        self.assertEqual(inf_q.versor, Quaternion(1, 0, -1, 0).versor)

        with self.assertRaises(ZeroDivisionError):
            zero.versor
    
    def test_angle(self):
        self.assertAlmostEqual(largest_q.angle, Quaternion(1, 1, 1, 1).angle, places=15)
        self.assertEqual(largest_q.vector.angle, 0.5*pi)
        self.assertAlmostEqual(q1.vector.angle, 0.5*pi, places=15)
        self.assertAlmostEqual(q2.vector.angle, 0.5*pi, places=15)
        self.assertEqual(q1.angle, q2.angle)
        self.assertEqual(zero.angle, 0.0)
        self.assertEqual(one.angle, 0.0)
        self.assertAlmostEqual(large_q.angle, random_q.angle, places=15)
        self.assertAlmostEqual(small_q.angle, random_q.angle, places=15)
        self.assertEqual(inf_q.angle, 0.25*pi)
        self.assertEqual((-inf_q).angle, 0.75*pi)


if __name__ == "__main__":
    unittest.main()
