#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes sure there are no obvious bugs in qmath functions.
"""
import sys
import unittest
from math import tau

from quaternions import qmath, Quaternion

large_float = sys.float_info.max
small_float = sys.float_info.min
TEN_TO_MINUS15 = 10e-15

# Main math function tests
q = Quaternion(1, -2, -3, 4)
exp_q = Quaternion(
    1.6939227236832994, 0.7895596245415588, 1.1843394368123383, -1.5791192490831176)
log_q = Quaternion(
    1.7005986908310777, -0.5151902926640851, -0.7727854389961277, 1.0303805853281702)
log10_q = Quaternion(
    0.7385606273598312, -0.2237443012341335, -0.33561645185120026, 0.447488602468267)
zero = Quaternion(0, 0, 0, 0)
one = Quaternion(1, 0, 0, 0)
ones = Quaternion(1, 1, 1, 1)


class TestQmathFunctions(unittest.TestCase):
    def test_exp(self):
        self.assertEqual(qmath.exp(1), qmath.e)
        self.assertEqual(qmath.exp(one), qmath.e)
        self.assertAlmostEqual(qmath.exp(q), exp_q, places=12)
    
    def test_log(self):
        self.assertEqual(qmath.log(1), 0)
        self.assertEqual(qmath.log(one), 0)
        self.assertAlmostEqual(qmath.log(q), log_q, places=12)
        self.assertAlmostEqual(qmath.log10(q), log10_q, places=12)
        
        with self.assertRaises(ValueError):
            qmath.log(zero)
    
    def test_sqrt(self):
        quaternions = [
            ones,
            Quaternion(1, -1, 1, 1), 
            Quaternion(1, 1, -1, 1), 
            Quaternion(1, 1, 1, -1),
        ]
        sizes = [2, 0.25, small_float, 0.5*large_float, large_float]
        for size in sizes:
            kwargs = dict(delta=size*TEN_TO_MINUS15)
            for q in quaternions:
                q = size*q
                self.assertAlmostEqual(qmath.sqrt(q)**2, q, **kwargs)
                self.assertAlmostEqual(qmath.sqrt(-q)**2, -q, **kwargs)
                self.assertAlmostEqual(qmath.sqrt(q.conjugate())**2, q.conjugate(), **kwargs)
                self.assertAlmostEqual(qmath.sqrt(-q.conjugate())**2, -q.conjugate(), **kwargs)
                self.assertAlmostEqual(qmath.sqrt(q.vector)**2, q.vector, **kwargs)
                self.assertAlmostEqual(qmath.sqrt(-q.vector)**2, -q.vector, **kwargs)
                self.assertAlmostEqual(qmath.sqrt(q.scalar)**2, q.real, **kwargs)
        
        self.assertAlmostEqual(
            qmath.sqrt(large_float*ones),
            large_float**0.5*qmath.sqrt(ones), 
            delta=large_float*TEN_TO_MINUS15
        )
        # Square root of negative real number.
        with self.assertRaises(ValueError):
            qmath.sqrt(-one)
        with self.assertRaises(ValueError):
            qmath.sqrt(-1.0)
    
    def test_rotate3d(self):
        p = (1,)
        self.assertEqual(qmath.rotate3d(p, 90, rounding=12), (0.0, 1.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 180, rounding=12), (-1.0, 0.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 270, rounding=12), (0.0, -1.0, 0.0))
        self.assertEqual(qmath.rotate3d(p, 360, rounding=12), (1.0, 0.0, 0.0))
    
    def test_cos(self):
        values = {
            tau/6: 0.5,
            tau/4: 0.0,
            tau/3: -0.5,
            tau/2: -1.0,
            3*tau/4: 0.0,
            tau: 1.0,
        }
        for i, o in values.items():
            kwargs = dict(delta=TEN_TO_MINUS15)
            self.assertAlmostEqual(qmath.cos(Quaternion(i)), Quaternion(o), **kwargs)
            self.assertAlmostEqual(qmath.cos(-Quaternion(i)), Quaternion(o), **kwargs)
        
        #TODO: make test cases for nonscalar quaternions.
    
    def test_sin(self):
        values = {
            tau/12: 0.5,
            tau/4: 1.0,
            5*tau/12: 0.5,
            tau/2: 0.0,
            3*tau/4: -1.0,
            tau: 0.0,
        }
        for i, o in values.items():
            kwargs = dict(delta=TEN_TO_MINUS15)
            self.assertAlmostEqual(qmath.sin(Quaternion(i)), Quaternion(o), **kwargs)
            self.assertAlmostEqual(qmath.sin(-Quaternion(i)), -Quaternion(o), **kwargs)
        
        #TODO: make test cases for nonscalar quaternions.
    
    #TODO: make test cases for newly added functions like cosh()


if __name__ == "__main__":
    unittest.main()
