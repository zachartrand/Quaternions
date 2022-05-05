Overview
========

The Quaternion Class
--------------------

The main aspect of Quaternions for Python is the Quaternion class.
The best way to use it is to import it like so:

>>> from quaternions import Quaternion

To create a quaternion, simply type

>>> Quaternion(a, b, c, d)

where a, b, c, and d correspond to a quaternion of the form :math:`a + bi + cj + dk`.
For example, creating the quaternion :math:`1 - 2i - 3j + 4k` looks like this in the
Python interpreter:


>>> q1 = Quaternion(1, -2, -3, 4)
>>> q1
Quaternion(1.0, -2.0, -3.0, 4.0)
>>> print(q1)
(1 - 2i - 3j + 4k)

Quaternions have mathematical functionality built in. Adding or multipling two
quaternions together uses the same syntax as ints and floats:


>>> q1, q2 = Quaternion(1, -2, -3, 4), Quaternion(1, 4, -3, -2)
>>> print(q1)
(1 - 2i - 3j + 4k)
>>> print(q2)
(1 + 4i - 3j - 2k)
>>> print(q1 + q2)
(2 + 2i - 6j + 2k)
>>> print(q1 - q2)
(-6i + 0j + 6k)
>>> print(q2 - q1)
(6i + 0j - 6k)
>>> print(q1 * q2)
(8 + 20i + 6j + 20k)
>>> print(q2 * q1)
(8 - 16i - 18j - 16k)
>>> print(q1/q2)
(-0.19999999999999996 - 0.8i - 0.4j - 0.4k)
>>> print(1/q2 * q1)
(-0.19999999999999996 + 0.4i + 0.4j + 0.8k)
>>> print(q2/q1)
(-0.19999999999999996 + 0.8i + 0.4j + 0.4k)


The qmath module
----------------

The qmath module functions similarly to Python's built-in :py:mod:`cmath` module
for complex numbers, allowing mathematical functions to be compatible with
quaternions. Here are a few examples:

>>> from quaternions import Quaternion, qmath
>>>
>>> q = Quaternion(1, -2, -3, 4)
>>> print(q)
(1 - 2i - 3j + 4k)
>>>
>>> print(qmath.exp(q))
(1.6939227236832994 + 0.7895596245415588i + 1.1843394368123383j - 1.5791192490831176k)
>>>
>>> print(qmath.log(q))
(1.7005986908310777 - 0.5151902926640851i - 0.7727854389961277j + 1.0303805853281702k)
>>>
>>> print(qmath.sqrt(q))
(1.7996146219471076 - 0.5556745248702426i - 0.833511787305364j + 1.1113490497404852k)
