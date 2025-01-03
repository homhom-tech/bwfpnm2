r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""

import scipy as _sp


def bulk_diffusion(physics,
                    phase,
                    network,
                    pc,
                    knudsen=True,
                    film_thickness='throat.film_thickness',
                    film_area='throat.film_area',
                    corner_area='throat.corner_area',
                    pore_length='throat.length', # or 'throat.porelengths'
                    pore_area='throat.area',
                    pore_diameter='throat.diameter',
                    pore_temperature='prop.temperature',
                    pore_waterdensity='prop.density',
                    pore_gasconstant='prop.gas_constant',
                    pore_pvsat='prop.Pvsat',
                    pore_permeability='prop.permeability',
                    eps=1e-12,
                    **kwargs):
    r"""
    Calculate the diffusive conductance of each element:
        kv = delta_v*Pv/(rho*Rv*T)/(1+Nk)*A/L
    where
        delta_v = Dab/(Rv*T)
        # vapour permeability in still air: (see diffusivity.py)
        Dab = 2.262/P*(T/273.15)**1.81  or  1.9e-10 kg/msPa at T=20C

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : VAPOUR Object
        The phase of interest

    Notes
    -----
    This function requires that all the necessary phase properties already be
    calculated.

    """
    _phases = physics._net._phases
    for _phase in _phases:
        if 'water' in _phase.name.split('_'):
            liquid = _phase
            break
    T = phase[pore_temperature][0]
    P = phase['prop.pressure'][0]
    rho = liquid[pore_waterdensity][0]
    Rv = phase[pore_gasconstant][0]
    Pvsat = phase[pore_pvsat][0]
    delta_v = phase[pore_permeability][0]

    # Find g for full throat
    if film_thickness.split('.')[0] == 'pore':
        # pore diffusive conductance
        conns = network['throat.conns']
        pore_length = 'throat.porelengths'
        pore_area = 'pore.area'
        pore_diameter = 'pore.diameter'
#        pore_shapefactor = 'pore.shapefactor'
#        pore_shapefactor_constant = 'pore.shapefactor_constant'
    else:
        conns = _sp.arange(network[pore_area].size)
#        pore_shapefactor = 'throat.shapefactor'
#        pore_shapefactor_constant = 'throat.shapefactor_constant'
    Apore = network[pore_area][conns]
    Lpore = network[pore_length]
    Lpore[Lpore <= 0] = eps

    rhoRT = rho*Rv*T
    rh = _sp.exp(pc/rhoRT)    # relative humidity
    Pv = Pvsat*rh             # vapour pressure
    kv = delta_v*Pv/rhoRT     # vapour permeability


    # Knudssen Diffusion
    if knudsen:
        rad = network[pore_diameter][conns]/2
        try:
            rad -= liquid[film_thickness][conns]
            rad[rad<=0] = eps
        except:
            pass
        Nk = knudsen_number(rad, T, P)
        kv = kv/(1+Nk)
#    # check if vapour conductance > liquid conductance [checked: never fulfilled]
#    G = network[pore_shapefactor][conns]
#    cG = network[pore_shapefactor_constant][conns]
#    mu = phase['prop.viscosity'][0]
#    Vapmax = delta_v*Pv*mu/rhoRT/rho
#    MatA = (1+Nk)*cG*G*Apore**2
    # Area constriction due to surface and corner adsorption
    try:
        Afilm = liquid[film_area][conns]
        Apore -= Afilm
    except:
        pass
    try:
        Acorner = liquid[corner_area][conns]
        Apore -= Acorner
    except:
        pass
#    MatA /= Apore
#    print('Vapmax & minMatA: {} vs {}'.format(Vapmax, MatA.min()))
#    try:
#        if pc > -10**5:
#            pass
#        tag = _sp.where(MatA < Vapmax)[0]
#        rtag = rad[tag]
#        Atag = network[pore_area][conns][tag]
#
#        if tag:
#            print('Diffusive conductance is higher than the hydraulic conductance\
#                   at pc {} at {} number {}'.format(pc, pore_area.split('.')[0], tag))
#    except:
#        pass

    gt = kv*Apore/Lpore          # conductance
    return gt


def knudsen_number(r, T, P, R=8.314, lm=None):
    r"""
    Quenard: lm = 100 nm = 1e-7
    Carmeliet: lm = 50 nm = 5e-8

    which one?
    applied: Drioli 2005
    """
    if lm is None:
        Mv, Ma = 18, 29 # g/mol
        Na = 6.0221e+23 # /mol
        dv, da = 0.27e-9, 0.37e-9 # m
        lm = R*T/(_sp.sqrt(1+Mv/Ma)*_sp.pi*Na*((dv+da)/2)**2*P)
    Kn = lm/(2*r)  # r: radius, lm: mean free path length
    return Kn


def tdiff_conduit(physics, throat_conductance='throat.diffusive_conductance',
                pore_conductance='throat.diffusive_conductance_pore',
                **kwargs):
    """
    Calculate vapour conduit \(1/2 pore - throat - 1/2 pore\) conductance.
    This function is used for calculating the singlephase permeability,
    e.g. called from singlephase_k.py.
    See Notebook on 2015-12-18.
    """
    gt = physics[throat_conductance]
    gp1 = physics[pore_conductance][:, 0]
    gp2 = physics[pore_conductance][:, 1]

    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return value
#def tbulk_diffusion_pore(physics,
#                         phase,
#                         network,
#                         pc,
#                         knudsen=True,
#                         throat_length='throat.porelengths',
#                         film_thickness='pore.surface_thickness_wp',
#                         pore_area='pore.area',
#                         pore_diameter='pore.diameter',
#                         pore_temperature='pore.temperature',
#                         pore_waterdensity='pore.water_density',
#                         pore_gasconstant='pore.gas_constant',
#                         pore_pvsat='pore.Pvsat',
#                         pore_permeability='pore.permeability',
#                         eps=1e-9,
#                         **kwargs):
#    r"""
#    Calculate the diffusive conductance of 2 pores connected to each throats. [Nt]
#    So the return gt dimension is (Nt, 2)
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    pc  : a constant
#
#    Notes
#    -----
#    This function requires that all the necessary phase properties already be
#    calculated.
#
#    """
#    T = phase[pore_temperature][0]
#    rho = phase[pore_waterdensity][0]
#    Rv = phase[pore_gasconstant][0]
#    Pvsat = phase[pore_pvsat][0]
#    delta_v = phase[pore_permeability][0]
#
#    # Find g for full throat
#    conns = network['throat.conns']
#    tarea = network[pore_area][conns]
#    tlen = network[throat_length]
#
#    #remove any non-positive lengths
#    tlen[tlen <= 0] = eps
#
#    rhoRT = rho*Rv*T
#    rh = _sp.exp(pc/rhoRT)    # relative humidity
#    Pv = Pvsat*rh             # vapour pressure
#    kv = delta_v*Pv/rhoRT     # vapour permeability
#
#    # Knudssen Diffusion
#    if knudsen:
#        rad = network[pore_diameter][conns]/2
#        try:
#            rad -= phase[film_thickness][conns]
#            rad[rad<=0] = eps
#            tarea = rad**2/network['pore.shapefactor'][conns]/4
#        except:
#            pass
#        Nk = knudsen_number(rad)
#        kv = kv/(1+Nk)
#
#    tarealen = tarea/tlen
#    gt = kv*tarealen          # conductance
#
#    return gt


#==============================================================================
# def pbulk_diffusion(physics,
#                     phase,
#                     network,
#                     pc,
#                     knudsen=True,
#                     pore_temperature='pore.temperature',
#                     pore_waterdensity='pore.water_density',
#                     pore_gasconstant='pore.gas_constant',
#                     pore_pvsat='pore.Pvsat',
#                     pore_permeability='pore.permeability',
#                     pore_area='pore.area',
#                     pore_diameter='pore.diameter',
#                     eps=1e-9,
#                     **kwargs):
#     r"""
#     Calculate the diffusive conductance of throats
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#
#     phase : OpenPNM Phase Object
#         The phase of interest
#
#     Notes
#     -----
#     This function requires that all the necessary phase properties already be
#     calculated.
#
#     """
#     T = phase[pore_temperature]
#     rho = phase[pore_waterdensity]
#     Rv = phase[pore_gasconstant]
#     Pvsat = phase[pore_pvsat]
#     Dab = phase[pore_permeability]
#
#     #Find g for full throat
#     parea = network[pore_area]
#     try:
#         plen = network[pore_diameter]/2
#     except:
#         plen = network['pore.length']
#
#     #remove any non-positive lengths
#     plen[plen <= 0] = eps
#
#     rhoRT = rho*Rv*T
#     rh = _sp.exp(pc/rhoRT)
#     Pv = Pvsat*rh
#     kv = Dab*Pv/rhoRT
#     gp = kv*(parea/plen)
#
#     # Knudssen Diffusion
#     if knudsen:
#         Nk = knudsen_number(plen)
#         gp = gp/(1+Nk)
#
#     return gp
#==============================================================================



#def tbulk_diffusion_eq(physics,
#                    phase,
#                    network,
#                    pc,
#                    knudsen=True,
#                    pore_temperature='pore.temperature',
#                    pore_waterdensity='pore.water_density',
#                    pore_gasconstant='pore.gas_constant',
#                    pore_pvsat='pore.Pvsat',
#                    pore_permeability='pore.permeability',
#                    throat_area='throat.area_eq',
#                    throat_length='throat.length',
#                    throat_radius='throat.radius_eq',
#                    film_thickness='throat.surface_thickness',
#                    eps=1e-9,
#                    **kwargs):
#    r"""
#    Calculate the diffusive conductance of throats
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    Notes
#    -----
#    This function requires that all the necessary phase properties already be
#    calculated.
#
#    """
#    T = phase[pore_temperature]
#    rhop = phase[pore_waterdensity]
#    Rvp = phase[pore_gasconstant]
#    Pvsat = phase[pore_pvsat]
#    delta_v = phase[pore_permeability]
#    #Interpolate pore phase property values to throats
#    T = phase.interpolate_data(data=T)
#    rhot = phase.interpolate_data(data=rhop)
#    Rvt = phase.interpolate_data(data=Rvp)
#    Pvsat = phase.interpolate_data(data=Pvsat)
#    delta_v = phase.interpolate_data(data=delta_v)   # vapour permeability
#
#    # Find g for full throat
#    tarea = network[throat_area]
#    tlen = network[throat_length]
#
#    #remove any non-positive lengths
#    tlen[tlen <= 0] = eps
#
#    rhoRT = rhot*Rvt*T
#    rh = _sp.exp(pc/rhoRT)    # relative humidity
#    Pv = Pvsat*rh             # vapour pressure
#    kv = delta_v*Pv/rhoRT     # vapour permeability
#    tarealen = tarea/tlen
#    gt = kv*tarealen          # conductance
#
#    # Knudssen Diffusion
#    if knudsen:
#        r = network[throat_radius] - phase[film_thickness]
#        Nk = knudsen_number(r)
#        gt = gt/(1+Nk)
#    return gt
#
#
#def pbulk_diffusion_eq(physics,
#                    phase,
#                    network,
#                    pc,
#                    knudsen=True,
#                    pore_temperature='pore.temperature',
#                    pore_waterdensity='pore.water_density',
#                    pore_gasconstant='pore.gas_constant',
#                    pore_pvsat='pore.Pvsat',
#                    pore_permeability='pore.permeability',
#                    pore_area='pore.area_eq',
#                    pore_radius='pore.radius_eq',
#                    pore_diameter='pore.diameter',
#                    film_thickness='pore.surface_thickness',
#                    eps=1e-9,
#                    **kwargs):
#    r"""
#    Calculate the diffusive conductance of throats
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    Notes
#    -----
#    This function requires that all the necessary phase properties already be
#    calculated.
#
#    """
#    T = phase[pore_temperature]
#    rho = phase[pore_waterdensity]
#    Rv = phase[pore_gasconstant]
#    Pvsat = phase[pore_pvsat]
#    Dab = phase[pore_permeability]
#
#    #Find g for full throat
#    parea = network[pore_area]
#    plen = network[pore_diameter]/2
#
#    #remove any non-positive lengths
#    plen[plen <= 0] = eps
#
#    rhoRT = rho*Rv*T
#    rh = _sp.exp(pc/rhoRT)
#    Pv = Pvsat*rh
#    kv = Dab*Pv/rhoRT
#    gp = kv*(parea/plen)
#
#    # Knudssen Diffusion
#    if knudsen:
#        # radius = pore radius - water film thickness
#        r = network[pore_radius] - phase[film_thickness]
#        Nk = knudsen_number(r)
#        gp = gp/(1+Nk)
#
#    return gp
#
#
#def tbulk_diffusion_pore_eq(physics,
#                         phase,
#                         network,
#                         pc,
#                         knudsen=True,
#                         pore_temperature='pore.temperature',
#                         pore_waterdensity='pore.water_density',
#                         pore_gasconstant='pore.gas_constant',
#                         pore_pvsat='pore.Pvsat',
#                         pore_permeability='pore.permeability',
#                         pore_area='pore.area_eq',
#                         throat_length='throat.porelengths',
#                         pore_radius='pore.radius_eq',
#                         film_thickness='pore.surface_thickness',
#                         eps=1e-9,
#                         **kwargs):
#    r"""
#    Calculate the diffusive conductance of pores connected to throats. [Nt]
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    pc  : a constant
#
#    Notes
#    -----
#    This function requires that all the necessary phase properties already be
#    calculated.
#
#    """
#    T = phase[pore_temperature][0]
#    rho = phase[pore_waterdensity][0]
#    Rv = phase[pore_gasconstant][0]
#    Pvsat = phase[pore_pvsat][0]
#    delta_v = phase[pore_permeability][0]
#
#    # Find g for full throat
#    conns = network['throat.conns']
#    tarea = network[pore_area][conns]
#    tlen = network[throat_length]
#    tfilm = phase[film_thickness][conns]
#
#    #remove any non-positive lengths
#    tlen[tlen <= 0] = eps
#
#    rhoRT = rho*Rv*T
#    rh = _sp.exp(pc/rhoRT)    # relative humidity
#    Pv = Pvsat*rh             # vapour pressure
#    kv = delta_v*Pv/rhoRT     # vapour permeability
#    tarealen = tarea/tlen
#    gt = kv*tarealen          # conductance
#
#    # Knudssen Diffusion
#    if knudsen:
#        r = network[pore_radius][conns] - tfilm
#        Nk = knudsen_number(r)
#        gt = gt/(1+Nk)              # shape (Nt, 2)
#    return gt


#def conduit_diffusion(physics, phase, network, pore_molar_density='pore.molar_density',
#                   pore_diffusivity='pore.diffusivity', pore_area='pore.area',
#                   pore_diameter='pore.diameter', throat_area='throat.area',
#                   throat_length='throat.length', throat_diameter='throat.diameter',
#                   shape_factor='throat.shape_factor', eps=1e-9, **kwargs):
#    r"""
#    Calculate the diffusive conductance of conduits in network, where a
#    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    Notes
#    -----
#    (1) This function requires that all the necessary phase properties already
#    be calculated.
#
#    (2) This function calculates the specified property for the *entire*
#    network then extracts the values for the appropriate throats at the end.
#
#    """
#    # Get Nt-by-2 list of pores connected to each throat
#    Ps = network['throat.conns']
#    # Get properties in every pore in the network
#    parea = network[pore_area]
#    pdia = network[pore_diameter]
#    # Get the properties of every throat
#    tdia = network[throat_diameter]
#    tarea = _sp.pi*(tdia/2)**2
#    tlen = network[throat_length]
#    # Interpolate pore phase property values to throats
#    cp = phase[pore_molar_density]
#    ct = phase.interpolate_data(data=cp)
#    DABp = phase[pore_diffusivity]
#    DABt = phase.interpolate_data(data=DABp)
#    if calc_pore_len:
#        lengths = misc.conduit_lengths(network, mode='centroid')
#        plen1 = lengths[:, 0]
#        plen2 = lengths[:, 2]
#    else:
#        plen1 = (0.5*pdia[Ps[:, 0]])
#        plen2 = (0.5*pdia[Ps[:, 1]])
#    # Remove any non-positive lengths
#    plen1[plen1 <= 0] = eps
#    plen2[plen2 <= 0] = eps
#    # Find g for half of pore 1
#    gp1 = ct*DABt*parea[Ps[:, 0]] / plen1
#    gp1[_sp.isnan(gp1)] = _sp.inf
#    gp1[~(gp1 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
#    # Find g for half of pore 2
#    gp2 = ct*DABt*parea[Ps[:, 1]] / plen2
#    gp2[_sp.isnan(gp2)] = _sp.inf
#    gp2[~(gp2 > 0)] = _sp.inf  # Set 0 conductance pores (boundaries) to inf
#    # Find g for full throat, remove any non-positive lengths
#    tlen[tlen <= 0] = eps
#    # Get shape factor
#    try:
#        sf = network[shape_factor]
#    except:
#        sf = _sp.ones(network.num_throats())
#    sf[_sp.isnan(sf)] = 1.0
#    gt = (1/sf)*ct*DABt*tarea/tlen
#    # Set 0 conductance pores (boundaries) to inf
#    gt[~(gt > 0)] = _sp.inf
#    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
#    value = value[phase.throats(physics.name)]
#    return value


