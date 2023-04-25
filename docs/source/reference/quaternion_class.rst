.. currentmodule:: quaternions

Quaternions
***********

The Quaternion class
====================
.. autoclass:: Quaternion

Methods
=======

Object-based methods
--------------------
.. automethod:: Quaternion.abs_components

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.from_complex

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.get_imag

.. automethod:: Quaternion.get_vector_components

Mathematical methods
--------------------
.. automethod:: Quaternion.conjugate

.. automethod:: Quaternion.inverse

.. automethod:: Quaternion.squared

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.cos_norm

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.sin_norm

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.log_norm

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.unit_quaternion

.. automethod:: Quaternion.unit_vector

Boolean methods
---------------
.. automethod:: Quaternion.is_complex

.. automethod:: Quaternion.is_scalar

.. automethod:: Quaternion.is_vector

.. automethod:: Quaternion.is_zero

   .. versionadded:: 2.0.0

.. automethod:: Quaternion.is_not_zero

   .. versionadded:: 2.0.0

Class methods
-------------
.. automethod:: Quaternion.from_angle

Properties
----------
.. autoproperty:: Quaternion.angle

.. property:: Quaternion.angle_in_radians

   Same as :py:attr:`Quaternion.angle`

.. autoproperty:: Quaternion.angle_in_degrees

.. autoproperty:: Quaternion.components

.. autoproperty:: Quaternion.norm

.. autoproperty:: Quaternion.scalar

.. autoproperty:: Quaternion.vector

.. autoproperty:: Quaternion.vector_norm

.. autoproperty:: Quaternion.versor
