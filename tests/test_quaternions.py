import unittest

from quaternions import Quaternion

q1 = Quaternion(1, -2, -3, 4)
q2 = Quaternion(1, 4, -3, -2)

q1_plus_q2 = Quaternion(2, 2, -6, 2)
q1_minus_q2 = Quaternion(0, -6, 0, 6)
q2_minus_q1 = Quaternion(0, 6, 0, -6)
q1_times_q2 = Quaternion(8.0, 20.0, 6.0, 20.0)
q2_times_q1 = Quaternion(8.0, -16.0, -18.0, -16.0)
q1_dividedby_q2 = Quaternion(-0.19999999999999996, -0.8, -0.4, -0.4)
q2_inverse_times_q1 = Quaternion(-0.19999999999999996, 0.4, 0.4, 0.8)
q2_dividedby_q1 = Quaternion(-0.19999999999999996, 0.8, 0.4, 0.4)
q1_inverse_times_q2 = Quaternion(-0.19999999999999996, -0.4, -0.4, -0.8)


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
        self.assertEqual(1/q2*q1, q2_inverse_times_q1)
        self.assertEqual(1/q1*q2, q1_inverse_times_q2)

    def test_division(self):
        self.assertEqual(q1/q2, q1_dividedby_q2)
        self.assertEqual(q2/q1, q2_dividedby_q1)

    def test_power(self):
        self.assertEqual(q1*q1, q1**2)
        self.assertEqual(q1.inverse(), q1**-1)
        self.assertEqual(q1, q1**1)
        self.assertEqual(q1**0, Quaternion(1))
        self.assertEqual(q1**-2, q1.inverse()*q1.inverse())


if __name__ == '__main__':
    unittest.main()
