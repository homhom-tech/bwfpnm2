# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 12:23:45 2015

@author: islah
"""
import scipy as sp


def constant(geometry, shapefactor, method='range', **kwargs):
    r'''
    Calculate constants based on given shape factors. These values are used
    in conductance calculation.

    Shape   G                     kG
    Circu    1/4/pi               0.5
    square    1/16                0.6
    triangle   (0, sqrt(3)/36]    0.5623

    Ref
    ---
    P. H. Valvatne, “Predictive pore-scale modelling of multiphase flow,” 2004.
    '''
    k = {'circular': 0.5, 'triangular': 0.6, 'square': 0.5623}
    G = {'circular': 1/4/sp.pi, 'square': 1/16, 'triangular': sp.sqrt(3)/36}

    if method == 'strict':
        circ = sp.isclose(shapefactor, G['circular'])
        squ = sp.isclose(shapefactor, G['square'])
        tri = ~(circ + squ)
    elif method == 'range':
        # based on Patzek & Silin 2001
        tri = shapefactor <= G['triangular']
        circ = shapefactor > G['square']
        squ = ~(tri + circ)

    k_G = sp.zeros_like(shapefactor)
    k_G[circ] = k['circular']
    k_G[squ] = k['square']
    k_G[tri] = k['triangular']
    return k_G

def shapefactor(geometry, radius, **kwargs):
    p = 2*sp.pi*radius
    area = sp.pi*radius**2
    G = area/p**2
    return G
