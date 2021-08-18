"""
Makes sure there are no obvious bugs in qmath functions.
"""
import unittest

from quaternions import qmath, Quaternion

# Main math function tests
q = Quaternion(1, -2, -3, 4)
exp_q = Quaternion(
    1.6939227236832994, 0.7895596245415588, 1.1843394368123383, -1.5791192490831176)
log_q = Quaternion(
    1.7005986908310777, -0.5151902926640851, -0.7727854389961277, 1.0303805853281702)
log10_q = Quaternion(
    0.7385606273598312, -0.2237443012341335, -0.33561645185120026, 0.447488602468267)


class TestQmathFunctions(unittest.TestCase):

    def test_exp(self):
        self.assertEqual(qmath.exp(1), qmath.e)
        self.assertEqual(qmath.exp(Quaternion(1)), qmath.e)
        self.assertEqual(qmath.exp(q), exp_q)

    def test_log(self):
        self.assertEqual(qmath.log(1), 0)
        self.assertEqual(qmath.log(Quaternion(1)), 0)
        self.assertEqual(qmath.log(q), log_q)
        self.assertEqual(qmath.log10(q), log10_q)

    def test_sqrt(self):
        self.assertEqual(qmath.sqrt(q), q**0.5)

    def test_rotate3d(self):
        p = (1,)
        self.assertEqual(qmath.rotate3d(p, 90, rounding=12), (0.0, 1.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 180, rounding=12), (-1.0, 0.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 270, rounding=12), (0.0, -1.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 360, rounding=12), (1.0, 0.0, 0.0))


if __name__ == '__main__':
    unittest.main()
