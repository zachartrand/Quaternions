import unittest

from math import pi

from quaternions import Quaternion

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


class TestQuaternions(unittest.TestCase):

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

    def test_division(self):
        self.assertAlmostEqual(q1/q2, q1_dividedby_q2, places=12)
        self.assertAlmostEqual(q2/q1, q2_dividedby_q1, places=12)

    def test_power(self):
        # Quaternion squared.
        self.assertEqual(q1*q1, q1**2)
        # Test pow(q, -1) == q.inverse()
        self.assertEqual(q1.inverse(), q1**-1)
        # q == pow(q, 1)
        self.assertEqual(q1, q1**1)
        # q**0 == 1
        self.assertEqual(q1**0, Quaternion(1))
        # pow(q, -2) == q.inverse()**2
        self.assertEqual(q1**-2, q1.inverse()*q1.inverse())

    def test_from_angle(self):
        # Test one case against expected values to 6 decimal places.
        self.assertAlmostEqual(
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False),
            Quaternion(2.500000, 1.157275, 2.314550, 3.471825), places=6)
        # Test the degrees flag as both True and False and make sure 
        # the result is the same.
        self.assertEqual(Quaternion.from_angle(60, [1, 2, 3], norm=5),
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False))
        # Test that an angle of zero produces a scalar quaternion.
        self.assertEqual(
            Quaternion.from_angle(0, [1, 2, 3], norm=5), Quaternion(5))
        # Test that an angle of zero produces a scalar quaternion 
        # even with an irrational norm.
        self.assertEqual(
            Quaternion.from_angle(0, [1, 2, 3], norm=pi),
            pi*Quaternion.from_angle(0, [1, 2, 3]))
        # Test that the norm parameter works properly, creating a 
        # quaternion of that norm.
        self.assertAlmostEqual(
            Quaternion.from_angle(pi/3, [1, 2, 3], norm=5, degrees=False).norm, 5.0, places=12)


if __name__ == "__main__":
    unittest.main()
