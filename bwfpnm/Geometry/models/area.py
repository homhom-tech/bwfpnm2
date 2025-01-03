# -*- coding: utf-8 -*-
r"""
===============================================================================
pore_area -- Models for cross-sectional area of a pore body
===============================================================================

"""
import scipy as _sp


def spherical(geometry,
              pore_diameter='pore.diameter',
              **kwargs):
    r"""
    Calculate cross-sectional area for a spherical pore
    """
    r = geometry[pore_diameter]/2
    value = _sp.pi*(r)**2
    return value


def cubic(geometry,
          pore_diameter='pore.diameter',
          **kwargs):
    r"""
    Calculate cross-sectional pore area for a cubic pore
    """
    diams = geometry[pore_diameter]
    value = (diams)**2
    return value


def shapefactor(geometry,
                shapefactor='pore.shapefactor',
                diameter='pore.diameter', **kwargs):

    r = geometry[diameter]/2
    G = geometry[shapefactor]

    A = r**2/G/4
    return A


def VL(geometry,
       volume='pore.volume',
       length='pore.length', **kwargs):

    if volume.split('.')[0] == 'pore':
        A = _VL_pore(geometry, volume, length, **kwargs)
    else:

        V = geometry[volume]
        L = geometry[length]

        A = V/L
    return A


def _VL_pore(geometry, volume, length, **kwargs):
    r''' Calculate the Area of pores '''
    net = geometry._net
    conns = net['throat.conns']
    pt = net.find_neighbor_throats(pores=net.pores(), flatten=False)

    V = geometry['pore.volume']
    Lt = geometry['throat.porelengths']

    At = V[conns]/Lt
    A = _sp.zeros_like(V)
    An = _sp.zeros_like(A)
    for i, con in enumerate(conns):
        A[con] += At[i]
        An[con] += 1
    A /= An
    return A
