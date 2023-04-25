"""
Quaternions for Python
======================

This package provides a class and mathematical operations for 
quaternions, an extension of the complex numbers.

This package's documentation assumes that you import the 
``Quaternion`` class like so:

>>> from quaternions import Quaternion

The other package is ``qmath``, a collection of mathematical 
functions that can be used on Quaternions, similar to Python's builtin
``cmath`` module for complex numbers. This package's documentation 
assumes that you import the package like so:

>>> from quaternions import qmath

To see all of the methods for the ``Quaternion`` class, type

>>> help(Quaternion)

To see all of the functions and constants within the ``qmath`` 
module, type

>>> help(qmath)

The rest of the functions and methods have thorough documentation on
what they do and how they are supposed to be used.

The main use case for Quaternions is 3d rotation, with the main 
function in this package being ``qmath.rotate3d``. However, the 
amount of Quaternion math programmed into this package is extensive, 
so feel free to use this to experiment and explore the nature of the
Quaternions!
"""

from .quaternions import *
