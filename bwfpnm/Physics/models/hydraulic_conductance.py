r"""
===============================================================================
Submodule -- hydraulic_conductance
===============================================================================

"""
import scipy as _sp
from bwfpnm.Geometry import models as gm


def hp_shapefactor(physics,
                   phase,
                   network,
                   pore_area='throat.area',
                   pore_length='throat.length',
                   pore_shapefactor='throat.shapefactor',
                   pore_shapefactor_constant='throat.shapefactor_constant',
                   pore_viscosity='prop.viscosity',
                   pore_density='prop.density',
                   eps=1e-12,
                   **kwargs):
    r"""
    Calculates the hydraulic conductivity of angular pore/throat
    using the shapefactore based Hagen-Poiseuille model.

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : WATER Object

    Ref
    -----
    [1] P. H. Valvatne, “Predictive pore-scale modelling of multiphase flow,” page 54, 2004.
    [2] M. Piri and M. Blunt, “Three-dimensional mixed-wet random pore-scale network modeling of two-and three-phase flow in porous media. I. Model description,” Phys. Rev. E, vol. 71, Feb. 2005.
    """
    mu = phase[pore_viscosity][0]
    rho = phase[pore_density][0]

    if pore_area.split('.')[0] == 'pore':
        # pore diffusive conductance
        conns = network['throat.conns']
        pore_length = 'throat.porelengths'
        pore_shapefactor = 'pore.shapefactor'
        pore_shapefactor_constant = 'pore.shapefactor_constant'
    else:
        conns = _sp.arange(network[pore_area].size)

    G = network[pore_shapefactor][conns]
    cG = network[pore_shapefactor_constant][conns]
    Apore = network[pore_area][conns]
    Lpore = network[pore_length]

    Lpore[Lpore <= 0] = eps

    gt = cG*rho*G*Apore**2/(mu*Lpore)
    return gt


def surface_cond(physics,
                 phase,
                 network,
                 base='shapefactor',
                 film_area='throat.film_area', #+case
                 corner_area='throat.corner_area', #+case
                 pore_occupancy='throat.occupancy',
                 pore_length='throat.length', #porelengths
                 pore_area='throat.area', # pore
                 pore_shapefactor='throat.shapefactor', #pore
                 pore_shapefactor_constant='throat.shapefactor_constant', #pore
                 pore_viscosity='prop.viscosity',
                 pore_density='prop.density',
                 eps=1e-12,
                 **kwargs):
    r"""
    Calculates the hydraulic conductivity of the surface flow
    """
    mu = phase[pore_viscosity][0]
    rho = phase[pore_density][0]

    if film_area.split('.')[0] == 'pore':
        # pore diffusive conductance
        conns = network['throat.conns']
        pore_length = 'throat.porelengths'
        pore_area = 'pore.area'
        pore_shapefactor = 'pore.shapefactor'
        pore_shapefactor_constant = 'pore.shapefactor_constant'
    else:
        conns = _sp.arange(network[pore_area].size)

    occupancy = ~_sp.bool8(phase[pore_occupancy])
    G = network[pore_shapefactor][conns]
    cG = network[pore_shapefactor_constant]*occupancy
    cG = cG[conns]
    Apore = network[pore_area][conns]
    Lpore = network[pore_length]
    Lpore[Lpore <= 0] = eps

    Asurf = _sp.zeros_like(Lpore)
    try:
        Asurf += phase[film_area][conns]
    except:
        pass
    try:
        Asurf += phase[corner_area][conns]
    except:
        pass

    P = _sp.sqrt(Apore/G)
    if base == 'shapefactor':
        # shape factor based
        G = Asurf/P**2
        cG = gm.shapefactor.constant(network._geometries[0], G)
        gt = cG*rho*G*Asurf**2/(mu*Lpore)
    else:
        # HP cylindrical based
#        Dh = 4*Asurf/P             # hydraulic diameter
#        rad = Dh/2 * occupancy
        rad = Asurf/P * occupancy[conns]   # hydraulic radius
        gt = rho*rad**2/(8*mu)*Asurf/Lpore
    return gt


def thp_conduit(physics, throat_conductance='throat.hydraulic_conductance',
                pore_conductance='throat.hydraulic_conductance_pore',
                **kwargs):
    """
    Calculate water conduit \(1/2 pore - throat - 1/2 pore\) conductance.
    This function is used for calculating the singlephase permeability,
    e.g. called from singlephase_k.py.
    See Notebook on 2015-12-18.
    """
    gt = physics[throat_conductance]
    gp1 = physics[pore_conductance][:, 0]
    gp2 = physics[pore_conductance][:, 1]

    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    return value


def hagen_poiseuille(physics,
                      phase,
                      network,
                      pore_area='throat.area',
                      pore_length='throat.length',
                      pore_diameter='throat.diameter',
                      pore_viscosity='prop.viscosity',
                      pore_density='prop.density',
                      eps=1e-12,
                      **kwargs):
    r"""
    Calculates the hydraulic conductivity of throat assuming cylindrical
    geometry using the Hagen-Poiseuille model

    Parameters
    ----------
    network : OpenPNM Network Object

    phase : WATER Object
    """
    mu = phase[pore_viscosity][0]
    rho = phase[pore_density][0]

    if pore_area.split('.')[0] == 'pore':
        conns = network['throat.conns']
        pore_length = 'throat.porelengths'
        pore_diameter = 'pore.diameter'
    else:
        conns = _sp.arange(network[pore_area].size)
    tdia = network[pore_diameter][conns]
    tarea = network[pore_area][conns]
    tlen = network[pore_length]
    tlen[tlen <= 0] = eps

    # gt below is equal to (pi*r^2)/L*(rho*r^2/8*mu)
    k = rho*(tdia/2)**2/(8*mu)
    gt = k*(tarea/tlen)

    return gt

#%% =========== DEPRECATED! uncommented due to old data compatibility ====
#==============================================================================
#
# def thp_shapefactor(physics,
#                    phase,
#                    network,
#                    pore_area='throat.area',
#                    pore_length='throat.length',
#                    pore_shapefactor='throat.shapefactor',
#                    pore_shapefactor_constant='throat.shapefactor_constant',
#                    pore_viscosity='prop.viscosity',
#                    pore_density='prop.density',
#                    eps=1e-12,
#                    **kwargs):
#     pass
#     return
# def thp_shapefactor_pore(physics,
#                          phase,
#                          network,
#                          pore_area='pore.area',
#                          throat_length='throat.porelengths',
#                          pore_shapefactor='pore.shapefactor',
#                          pore_shapefactor_constant='pore.shapefactor_constant',
#                          pore_viscosity='prop.viscosity',
#                          pore_density='prop.density',
#                          eps=1e-12,
#                          **kwargs):
#     r"""
#     Calculates the hydraulic conductivity of pores connected to throat
#     assuming pore lengths based on Dong; Oren & Bakke, which are supplied in
#     the topological data.
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#
#     phase : OpenPNM Phase Object
#
#     Ref
#     -----
#     [1] P. H. Valvatne, “Predictive pore-scale modelling of multiphase flow,” page 54, 2004.
#     [1] M. Piri and M. Blunt, “Three-dimensional mixed-wet random pore-scale network modeling of two-and three-phase flow in porous media. I. Model description,” Phys. Rev. E, vol. 71, Feb. 2005.
#     """
#     mu = phase[pore_viscosity][0]
#     rho = phase[pore_density][0]
#
#     conns = network['throat.conns']
#     tshapefactor = network[pore_shapefactor][conns]
#     k_G = network[pore_shapefactor_constant][conns]
#     tarea = network[pore_area][conns]
#     tlen = network[throat_length]
#
#     # remove any non-positive lengths
#     tlen[tlen <= 0] = eps
#
#     gt = k_G*tshapefactor*tarea**2/mu
#     gt = rho*gt/tlen
#
#     return gt
#
#
#
#
#
# def phagen_poiseuille(physics,
#                       phase,
#                       network,
#                       pore_diameter='pore.diameter',
#                       pore_viscosity='prop.viscosity',
#                       pore_density='prop.density',
#                       pore_area='pore.area',
#                       eps=1e-9,
#                       **kwargs):
#     r"""
#     Calculates the hydraulic conductivity of pore assuming cylindrical
#     geometry using the Hagen-Poiseuille model
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#
#     phase : OpenPNM Phase Object
#     """
#     mup = phase[pore_viscosity]
#     rhop = phase[pore_density]
#     pdia = network[pore_diameter]
#     plen = 0.5*pdia             # network[pore_length]
#     parea = network[pore_area]
#     # remove any non-positive lengths
#     plen[plen <= 0] = eps
#
#     # gp below is equal to (pi*r^2)/L*(rho*r^2/8*mu)
#     k = rhop*(pdia/2)**2/(8*mup)
#     gp = k*(parea/plen)
#     return gp
#
#
#
#
#
# def php_shapefactor(physics,
#                     phase,
#                     network,
#                     pore_diameter='pore.diameter',
#                     pore_viscosity='prop.viscosity',
#                     pore_density='prop.density',
#                     pore_area='pore.area',
#                     pore_shapefactor='pore.shapefactor',
#                     pore_shapefactor_constant='pore.shapefactor_constant',
#                     eps=1e-9,
#                     **kwargs):
#     r"""
#     Calculates the hydraulic conductivity of pore assuming pore length = radius
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#
#     phase : OpenPNM Phase Object
#     """
#     pshapefactor = network[pore_shapefactor]
#     k_G = network[pore_shapefactor_constant]
#     mup = phase[pore_viscosity]
#     rhop = phase[pore_density]
#     parea = network[pore_area]
#
#     # should we use throat's porelength? YES! see below 'thp_shapefactor_pore'
#     plen = 0.5*network[pore_diameter]
#
#     # remove any non-positive lengths
#     plen[plen <= 0] = eps
#
#     gp = k_G*pshapefactor*parea**2/mup
#     gp = rhop*gp/plen
#
#     return gp
#
#==============================================================================
