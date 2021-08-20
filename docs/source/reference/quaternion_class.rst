.. currentmodule:: quaternions

Quaternions
***********

The Quaternion class
====================
.. autoclass:: Quaternion

Methods
=======

Normal methods
--------------
.. automethod:: Quaternion.conjugate

.. automethod:: Quaternion.get_imag

.. automethod:: Quaternion.get_vector_components

.. automethod:: Quaternion.inverse

.. automethod:: Quaternion.unit_quaternion

.. automethod:: Quaternion.unit_vector

Boolean methods
---------------
.. automethod:: Quaternion.is_complex

.. automethod:: Quaternion.is_scalar

.. automethod:: Quaternion.is_vector

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
