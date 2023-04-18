"""
Module containing miscellaneous helper functions and definitions.
"""

from math import (
    copysign as _copysign,
    cos as _cos,
    sin as _sin,
    radians as _radians,
    inf as _inf,
    nan as _nan,
    isfinite as _isfinite,
    isnan as _isnan,
)
from typing import Iterable


def makeListLen3(i: Iterable[float]) -> list:
    """Makes sure points and axes of rotation have 3 coordinates."""
    if not isinstance(i, list):
        i = list(i)

    if len(i) > 3:
        i = i[:3]
    elif len(i) < 3:
        while len(i) < 3:
            i.append(0.0)
    return i


def deg_cos(x: float) -> float:
    """Return cos(x) where x is in degrees."""
    EXACTS = {
        0.0: 1.0, 
        60.0: 0.5, 
        90.0: 0.0, 
        120.0: -0.5, 
        180.0: -1.0,
        240.0: -0.5,
        270.0: 0.0,
        300.0: 0.5,
    }
    if x < 0.0:
        x = abs(x)
    
    reduced_x = x % 360.0
    
    return EXACTS.get(reduced_x, _cos(_radians(x)))


def deg_sin(x: float) -> float:
    """Return sin(x) where x is in degrees."""
    EXACTS = {
        0.0: 0.0, 
        30.0: 0.5, 
        90.0: 1.0, 
        150.0: 0.5, 
        180.0: 0.0,
        210.0: -0.5,
        270.0: -1.0,
        330.0: -0.5,
    }
    negative = 1.0
    if x < 0.0:
        negative = -1.0
    
    x = abs(x)
    reduced_x = x % 360.0
    
    return negative*EXACTS.get(reduced_x, _sin(_radians(x)))
