# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 14:50:20 2015

@author: islah
"""
import scipy as _sp


def stat_thickness(phase, pc,
                   K1=3.85, K2=-1.89,
                   Rv=462,
                   trapping=False,
                   trap_pc='pore.imbibition_trapped_pc',
                   film_thickness='pore.film_thickness',
                   pore_occupancy='pore.occupancy_wp',
                   pore_temperature='prop.temperature',
                   pore_waterdensity='prop.density',
                   pore_gasconstant='prop.gas_constant',
                   **kwargs):
    r"""
    Calculate the statistical thickness of water film on pore wall at specified
    pc/relative humidity based on Bradley's equation.

    The thickness will be negative below around 0.07% relative humidity, and
    goes to infinity for 100% RH.

    Parameters
    ----------
    phase : WATER Object
        The phase of interest

    K1, K2 : material parameters, for silicate materials: K1=3.85, K2=-1.89

    Ref
    -----
    [1] R. Badmann, N. Stockhausen, and M. J. Setzer, “The statistical thickness
    and the chemical potential of adsorbed water films,” J. Colloid Interface
    Sci., vol. 82, no. 2, pp. 534–542, 1981.

    """
    T = phase[pore_temperature][0]
    rho = phase[pore_waterdensity][0]
    occupancy = ~_sp.bool8(phase[pore_occupancy])
    rad = phase._net['pore.diameter']/2

    if pore_occupancy.split('.')[0] == 'throat':
        rad = phase._net['throat.diameter']/2

    if trapping:
        p_trap = pc > phase[trap_pc]        # entrapped pores
    else:
        p_trap = _sp.zeros_like(rad, dtype=bool)
    p_actv = _sp.bool8(occupancy*1 - p_trap*1)         # active dry pores
#    p_actv = occupancy - p_trap        # active dry pores
    try:
        # the old corner area for the entrapped pores
        t_old = phase[film_thickness]
    except:
        t_old = _sp.zeros_like(rad)

    # statistical thickness
#    rh = _sp.exp(pc/(rho*Rv*T))
    t_new = (K1 + K2*_sp.log(-pc/(rho*Rv*T)))*1e-10  # a scalar

    if t_new<0.0:
        t_new = 0.0
    else:
        # check whether the thickness > radius
        tbig = t_new>rad
        t_new = t_new*occupancy
        if _sp.any(tbig):
            print('Warning: film thickness > pore radius at some pores!\n\
                  lPc: {}'.format(_sp.log10(-pc)))
            tag = _sp.where(tbig)[0]
            rat = _sp.sum(t_new[tag])/_sp.sum(rad[tag])
            print('Ratio (film/pore): {}'.format(rat))
        t_new[tbig] = rad[tbig]

    t = t_old * p_trap + t_new * p_actv
    return t


def area(phase,
         pore_occupancy='pore.occupancy_wp',
         film_thickness='pore.film_thickness_wp',
         pore_diameter='pore.diameter',
         pore_area='pore.area',
         **kwargs):
    r"""
    Calculate the area of the water films on pore wall at specified
    pc/relative humidity.

    Parameters
    ----------
    network : WATER Object

    phase : OpenPNM Phase Object
        The phase of interest

    """
    occupancy = ~_sp.bool8(phase[pore_occupancy])
    t_film = phase[film_thickness]
    if pore_occupancy.split('.')[0] == 'throat':
        pore_diameter = 'throat.diameter'
        pore_area = 'throat.area'
    r_pore = phase._net[pore_diameter]/2
    Apore = phase._net[pore_area]

    #%% surface film area
    # concruency: A_film = A_pore*(1 - (r_pore - t_film)**2/r_pore**2)
    Afilm = Apore*(1 - (r_pore - t_film)**2/r_pore**2)*occupancy

    Afilm[Afilm < 0] = 0.0
    Abig = Afilm > Apore
    if _sp.any(Abig):
        Afilm[Abig] = Apore[Abig]
    return Afilm


def volume(phase,
           pore_occupancy='pore.occupancy_wp',
           film_area='pore.film_area_wp',
           pore_volume='pore.volume',
           pore_area='pore.area',
           **kwargs):
    r"""
    Calculate the volume of the water films of pore wall at specified
    pc/relative humidity.

    Parameters
    ----------
    phase : WATER Object
        The phase of interest

    """
#    occupancy = phase[pore_occupancy]
    A_film = phase[film_area]
    if pore_occupancy.split('.')[0] == 'throat':
        pore_area = 'throat.area'
        pore_volume = 'throat.volume'
    V_pore = phase._net[pore_volume]
    A_pore = phase._net[pore_area]

    #%% surface film volume
    Vfilm = A_film/A_pore * V_pore  # * occupancy is implicit in A_film

    Vfilm[Vfilm < 0] = 0.0
    Vbig = Vfilm > V_pore
    if _sp.any(Vbig):
        Vfilm[Vbig] = V_pore[Vbig]
    return Vfilm

#def volume_eq(phase,
#              shapefactor='pore.shapefactor',
#              radius='pore.radius_eq',
#              length='pore.length',
#              volume='pore.volume',
#              occupancy='pore.occupancy_wp',
#              film_thickness='pore.surface_thickness_wp',
#              **kwargs):
#    r"""
#    Calculate the volume of the water films in pores/throats at specified
#    pc/relative humidity.
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    """
#    G = phase._net[shapefactor]
#    R = phase._net[radius]
#    V = phase._net[volume]
#    t = phase[film_thickness]
#    try:
#        L = phase._net[length]
#    except:
#        L = phase._net['pore.diameter']     # see when calculating area_eq
#    occupancy = phase[occupancy]
#
#    # surface film volume
##==============================================================================
##     r = R - t
##     Avapour = r**2/4/G
##     Vvapour = Avapour*L
##     Vfilm = (V - Vvapour)*occupancy
##==============================================================================
#
#    Vfilm = V*(t/R*(2-t/R))*occupancy
#    Vfilm[Vfilm < 0] = 0.0
#    return Vfilm
#
#
#def php_shapefactor(physics,
#                    phase,
#                    network,
#                    pore_diameter='pore.diameter',
#                    pore_viscosity='pore.viscosity',
#                    pore_density='pore.density',
#                    pore_area='pore.area',
#                    pore_shapefactor='pore.shapefactor',
#                    pore_shapefactor_constant='pore.shapefactor_constant',
#                    eps=1e-9,
#                    **kwargs):
#    r"""
#    Calculates the hydraulic conductivity of pore assuming pore length = radius
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#    """
#    pshapefactor = network[pore_shapefactor]
#    k_G = network[pore_shapefactor_constant]
#    mup = phase[pore_viscosity]
#    rhop = phase[pore_density]
#    parea = network[pore_area]
#
#    # should we use throat's porelength? YES! see below 'thp_shapefactor_pore'
#    plen = 0.5*network[pore_diameter]
#
#    # remove any non-positive lengths
#    plen[plen <= 0] = eps
#
#    gp = k_G*pshapefactor*parea**2/mup
#    gp = rhop*gp/plen
#
#    return gp
#
#    # surface film conductance
##==============================================================================
##     def conductance_water(r, A, L, G, c_G, rho=1000, mu=1e-3):
##     # cylindrical pore
##     glo = rho*r**2/(8*mu) * A/L
##
##     # arbitrary-shape pore: the shape factor is not relevant anymore since
##     # the area changes although the perimeter is still ok.
##     gl = c_G*rho*G*A/mu * A/L
##==============================================================================
