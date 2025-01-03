# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 04:30:18 2015

@author: islah

Calculate the size (radius) of the equivalent circular tube
conserving the real volume. This is necessary for moisture content calculation
from the surface adsorption (water films).
"""

import scipy as _sp


def radius(geometry,
           shape='circular',
           volume='pore.volume',
           length='pore.length',
           shapefactor='pore.shapefactor',
           **kwargs):
    r"""
    Calculate the equivalent radius of circular/arbitrary tube for conserving
    the same real volume of arbitrary shapes.
    """
    try:
        L = geometry[length]
    except:
        L = geometry['pore.diameter']

    A = geometry[volume]/L

    if shape == 'circular':
        r = _sp.sqrt(A/_sp.pi)      # assuming circular cross-section
    else:
        G = geometry[shapefactor]
        r = 2*_sp.sqrt(G*A)         # assuming arbitrary shapes, based on G
    return r


def area(geometry,
         volume='pore.volume',
         length='pore.length',
         **kwargs):
    r"""
    Calculate the equivalent area of arbitrary tube for conserving
    the same real volume of arbitrary shapes.
    """
    try:
        L = geometry[length]
    except:
        L = geometry['pore.diameter']

    A = geometry[volume]/L
    return A
