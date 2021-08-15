# Quaternions

Class and mathematical functions for quaternion numbers.

## Installation

### Python

This is a Python 3 module.  If you don't have Python installed, get the latest
version [here](https://www.python.org/downloads/).

### The Quaternions module

Install with pip:
```
pip install quaternions-for-python
```

If you want to build from source, you can clone the repository with the following
terminal command:
```
git clone https://github.com/zachartrand/Quaternions.git
```

## How to use

### Using the quaternions module

The quaternions module is designed to be imported to use quaternion numbers
just like complex numbers in Python. The rest of this webpage assumes you
import the class like this:

```python
>>> from quaternions import Quaternion
```

To create a quaternion, simply type
```python
>>> Quaternion(a, b, c, d)
```
where a, b, c, and d correspond to a quaternion of the form `a + bi + cj + dk`.
For example, creating the quaternion `1 - 2i - 3j + 4k` looks like this in the
Python interpreter:

```python
>>> q1 = Quaternion(1, -2, -3, 4)
>>> q1
Quaternion(1.0, -2.0, -3.0, 4.0)
>>> print(q1)
(1 - 2i - 3j + 4k)
```

Quaternions have mathematical functionality built in. Adding or multipling two
quaternions together uses the same syntax as ints and floats:

```python
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
```

Check the documentation for other useful methods of the Quaternion class.

### Using the qmath module
The qmath module contains some functions that are compatible with quaternions,
similarly to how the cmath module works. These include the exponential function,
the natural logarithm, and the pow function. It also includes a function,
rotate3d, that takes an iterable of coordinates and rotates them a given angle
around a given axis (the z-axis by default). Here is an example rotating the
point (1, 0, 0) around the z-axis:
```python
>>> from quaternions import qmath
>>>
>>> p = (1, 0, 0)
>>>
>>> p = qmath.rotate3d(p, 90); print(p)
(0.0, 1.0, 0.0)
>>> p = qmath.rotate3d(p, 90); print(p)
(-1.0, 0.0, 0.0)
>>> p = qmath.rotate3d(p, 90); print(p)
(0.0, -1.0, 0.0)
>>> p = qmath.rotate3d(p, 90); print(p)
(1.0, 0.0, 0.0)
```
