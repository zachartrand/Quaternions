The qmath module
================

qmath
-----
.. automodule:: quaternions.qmath
  :no-members:

.. currentmodule:: quaternions.qmath

Quaternion mathematical functions
---------------------------------
.. autofunction:: exp

.. autofunction:: log

.. autofunction:: log10

.. autofunction:: sqrt

Rotation functions
------------------
.. autofunction:: rotate3d

.. autofunction:: rotate_Euler

   .. versionadded:: 1.1.0

Other functions
---------------
.. autofunction:: cross_product

   .. versionadded:: 1.1.0

.. autofunction:: dot_product

   .. versionadded:: 1.1.0

Constants
---------
.. data:: pi

   The mathematical constant *π*, as a float.

.. data:: tau

   The mathematical constant *τ*, as a float.

.. data:: e

   The mathematical constant *e*, as a float.

.. data:: inf

   Floating-point positive infinity. Equivalent to ``float('inf')``.

.. data:: infi

   Quaternion with positive infinity i part and zero for all the other parts.
   Equivalent to ``Quaternion(0, float('inf'), 0, 0)``.

.. data:: infj

  Quaternion with positive infinity j part and zero for all the other parts.
  Equivalent to ``Quaternion(0, 0, float('inf'), 0)``.

.. data:: infk

  Quaternion with positive infinity k part and zero for all the other parts.
  Equivalent to ``Quaternion(0, 0, 0, float('inf'))``.

.. data:: nan

   A floating-point “not a number” (NaN) value. Equivalent to ``float('nan')``.

.. data:: nani

      Quaternion with NaN i part and zero for all the other parts.
      Equivalent to ``Quaternion(0, float('nan'), 0, 0)``.

.. data:: nanj

   Quaternion with NaN j part and zero for all the other parts.
   Equivalent to ``Quaternion(0, 0, float('nan'), 0)``.

.. data:: nank

   Quaternion with NaN k part and zero for all the other parts.
   Equivalent to ``Quaternion(0, 0, 0, float('nan'))``.
