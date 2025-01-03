#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 16:51:00 2017

@author: islah
"""
import scipy as _sp


def area(phase, pc,
         trapping=False,
         trap_pc='pore.imbibition_trapped_pc',
         film_thickness='pore.film_thickness',
         film_area='pore.film_area',
#         film_vol='pore.film_volume',
         corner_area='pore.corner_area',
         pore_occupancy='pore.occupancy_wp',
         pore_area='pore.area',
         pore_shapefactor='pore.shapefactor',
         pore_surfacetension='prop.surface_tension',
         **kwargs):
    r"""
    Calculate the area of corner water in dry angular pores (square and
    triangular cross sections).

    Parameters
    ----------

    Ref
    -----
    [1] G. Mason and N. R. Morrow, “Capillary behavior of a perfectly
    wetting liquid in irregular triangular tubes,” J. Colloid Interface Sci.,
    vol. 141, no. 1, pp. 262–274, Jan. 1991.
    """
    sigma = phase[pore_surfacetension][0]
    element = pore_occupancy.split('.')[0]
    if element == 'throat':
        pore_area = 'throat.area'
        pore_shapefactor = 'throat.shapefactor'
    G = phase._net[pore_shapefactor]
    A = phase._net[pore_area]
    P = _sp.sqrt(A/G)

    Rx = phase._net[element+'.diameter']/2 - phase[film_thickness]
    Rpc = -(1+2*_sp.sqrt(_sp.pi*G))*sigma*_sp.cos(0)/pc
    ind = Rpc > Rx    # pore indices where the r_corner reaches the max rad

    try:
        # the old corner area for the entrapped pores
        Ac_old = phase[corner_area]
    except:
        Ac_old = _sp.zeros_like(A)

    # pore shapes: corner ad only occurs in angular shapes
    squ = phase._net[element+'.square']
    tri = phase._net[element+'.triangular']

    # adjust pore size due to water film
    try:
        Ax = A - phase[film_area]   # constricted area = pore area - film area
        P *= _sp.sqrt(Ax/A)         # Pf = k*P, Ax = A-At = k^2 * A
    except:
        Ax = A
    # =====================
    # small pc ==> small r ==> small corner area
    # critical pc ==> crit r_d = r_inscribed ==> max corner area
    # pc > crit pc ==> condensation ==> Sat = 1
    # So, the following Sw equation only valid for r <= r_d = r_inscribed
    # This condition is automatically fulfilled in wetting & drying cases,
    # but not in imbibition with air entrapment! Hence, the corner adsorption
    # must not occurs once the pores are entrapped.
    # =====================
    # dry pores (hence we Ac[wet pores] = 0)
    occupancy = ~_sp.bool8(phase[pore_occupancy])   # dry pores
    # entrapped dry pores
    if trapping:
        p_trap = pc > phase[trap_pc]        # entrapped pores
    else:
        p_trap = _sp.zeros_like(G, dtype=bool)
    p_actv = _sp.bool8(occupancy*1 - p_trap*1)         # active dry pores
#    p_actv = occupancy - p_trap         # active dry pores

    # saturation of corner adsorption
    Sw = _sp.zeros_like(G)      # Sw[circ] = 0
    Sw[tri] = (Rpc[tri]/P[tri])**2 *(1/G[tri])*(0.25/G[tri] - _sp.pi)
    Sw[squ] = ((4 - _sp.pi)*Rpc[squ]**2)/Ax[squ]

    sw1 = Sw > 1
    Sw[ind] = 1 # w.r.t. the constricted area Ax (see the next command line)
    Sw[sw1] = 1

    # area of corner adsorption
    Ac_new = Sw * Ax # * occupancy     # here Ac[wet_pores] = 0
    Ac = Ac_old * p_trap + Ac_new * p_actv
    if _sp.any(Ac > Ax):
        print('Warning: corner area > pore area at some pores!\n\
              lPc: {}'.format(_sp.log10(-pc)))
        tag = _sp.where(Ac > Ax)[0]
        rat = _sp.sum(Ac[tag])/_sp.sum(Ax[tag])
        print('Ratio (corner/pore): {}'.format(rat))
    return Ac

def volume(phase,
           pore_occupancy='pore.occupancy_wp',
           corner_area='pore.corner_area',
           pore_area='pore.area',
           pore_volume='pore.volume',
           **kwargs):
    r"""
    Calculate the volume of corner water in dry angular pores (square and
    triangular cross sections).

    Parameters
    ----------

    Ref
    -----
    [1] G. Mason and N. R. Morrow, “Capillary behavior of a perfectly
    wetting liquid in irregular triangular tubes,” J. Colloid Interface Sci.,
    vol. 141, no. 1, pp. 262–274, Jan. 1991.

    """
#    occupancy = phase[pore_occupancy]
    if pore_occupancy.split('.')[0] == 'throat':
        pore_area = 'throat.area'
        pore_volume = 'throat.volume'
    V = phase._net[pore_volume]
    A = phase._net[pore_area]
    Ac = phase[corner_area]

    # corner volume
    Vc = Ac/A * V   # *~_sp.bool8(occupancy) is implicit in Ac
    return Vc
