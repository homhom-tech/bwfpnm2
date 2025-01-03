r"""
===============================================================================
Submodule -- capillary_pressure
===============================================================================

"""

import scipy as _sp

def washburn(physics,
             phase,
             network,
             pore_diameter='pore.diameter', #throat
             pore_surface_tension='prop.surface_tension',
             pore_contact_angle='prop.contact_angle',
             method='r-based', # or G-based
             ndec=3,
             pc_eps=-1e+12,
             **kwargs):
    r"""
    Computes the capillary entry pressure assuming the throat is a cylindrical tube.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network on which to apply the calculation
    phase : OpenPNM Phase Object
        Phase object for the invading phases

    Notes
    -----
    The Washburn equation is:

    .. math::
        P_c = -\frac{2\sigma(cos(\theta))}{r}

    This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading phases in most materials.

    """
    shapefactor = 'pore.shapefactor'
    sigma = phase[pore_surface_tension][0]
    theta = phase[pore_contact_angle][0]
    rad = network[pore_diameter]/2
    if pore_diameter.split('.')[0] == 'throat':
        poreloc = phase.throats(physics.name)   # assigned pores
        shapefactor = 'throat.shapefactor'
    else:
        poreloc = phase.pores(physics.name)

    if method.lower() == 'r-based':
        cc = 2     # for all pores
    elif method.lower() == 'g-based':
        G = network[shapefactor]
        cc = 1 + 2*_sp.sqrt(_sp.pi*G)

    value = -cc*sigma*_sp.cos(_sp.radians(theta))/rad

    # set pc=pc_eps when pc=-inf due to zero diameter
    x = _sp.amin(_sp.r_[value[value!=_sp.inf], pc_eps])
    value[value==-_sp.inf] = x
    value = value[poreloc]                         # only the assigned throats
    return _sp.around(value, ndec)

#%% =========== DEPRECATED! uncommented due to old data compatibility ====
#==============================================================================
#
# def twashburn(physics,
#               phase,
#               network,
#               pore_surface_tension='prop.surface_tension',
#               pore_contact_angle='prop.contact_angle',
#               pore_diameter='pore.diameter',
#               ndec=3,
#               **kwargs):
#     return
#
# def pwashburn(physics,
#               phase,
#               network,
#               pore_surface_tension='prop.surface_tension',
#               pore_contact_angle='prop.contact_angle',
#               pore_diameter='pore.diameter',
#               ndec=3,
#               **kwargs):
#     r"""
#     Computes the capillary entry pressure assuming the pore is a cylindrical tube.
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#         The network on which to apply the calculation
#     phase : OpenPNM Phase Object
#         Phase object for the invading phases
#
#     Notes
#     -----
#     The Washburn equation is:
#
#     .. math::
#         P_c = -\frac{2\sigma(cos(\theta))}{r}
#
#     This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading phases in most materials.
#
#     """
#     pores = phase.pores(physics.name)
#     sigma = phase[pore_surface_tension]
#     theta = phase[pore_contact_angle]
#     r = network[pore_diameter]/2
#
#     value = -2*sigma*_sp.cos(_sp.radians(theta))/r
#     # set pc=-1e+10 when pc=-inf due to zero diameter
#     x = _sp.amin(_sp.r_[value[value!=-_sp.inf], -1e+9])
#     value[value==-_sp.inf] = x
#     value = value[pores]
#     return _sp.around(value, ndec)
#
#
# def twashburn_eq(physics,
#                  phase,
#                  network,
#                  pore_surface_tension='prop.surface_tension',
#                  pore_contact_angle='prop.contact_angle',
#                  throat_radius='throat.radius_eq',
#                  ndec=3,
#                  **kwargs):
#     r"""
#     Computes the capillary entry pressure assuming the throat is a cylindrical tube.
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#         The network on which to apply the calculation
#     phase : OpenPNM Phase Object
#         Phase object for the invading phases
#
#     Notes
#     -----
#     The Washburn equation is:
#
#     .. math::
#         P_c = -\frac{2\sigma(cos(\theta))}{r}
#
#     This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading phases in most materials.
#
#     """
#     throats = phase.throats(physics.name)   # assigned throats
#     sigma = phase[pore_surface_tension]
# #    sigma = phase.interpolate_data(data=sigma)      # for all throats
#     theta = phase[pore_contact_angle]
# #    theta = phase.interpolate_data(data=theta)      # for all throats
#     r = network[throat_radius]
#
#     value = -2*sigma*_sp.cos(_sp.radians(theta))/r      # for all throats
#     # set pc=-1e+10 when pc=-inf due to zero diameter
#     x = _sp.amin(_sp.r_[value[value!=_sp.inf], -1e+9])
#     value[value==-_sp.inf] = x
#     value = value[throats]                         # only the assigned throats
#     return _sp.around(value, ndec)
#
#
# def pwashburn_eq(physics,
#                  phase,
#                  network,
#                  pore_surface_tension='prop.surface_tension',
#                  pore_contact_angle='prop.contact_angle',
#                  pore_radius='pore.radius_eq',
#                  ndec=3,
#                  **kwargs):
#     r"""
#     Computes the capillary entry pressure assuming the pore is a cylindrical tube.
#
#     Parameters
#     ----------
#     network : OpenPNM Network Object
#         The network on which to apply the calculation
#     phase : OpenPNM Phase Object
#         Phase object for the invading phases
#
#     Notes
#     -----
#     The Washburn equation is:
#
#     .. math::
#         P_c = -\frac{2\sigma(cos(\theta))}{r}
#
#     This is the most basic approach to calculating entry pressure and is suitable for highly non-wetting invading phases in most materials.
#
#     """
#     pores = phase.pores(physics.name)
#     sigma = phase[pore_surface_tension]
#     theta = phase[pore_contact_angle]
#     r = network[pore_radius]
#
#     value = -2*sigma*_sp.cos(_sp.radians(theta))/r
#     # set pc=-1e+10 when pc=-inf due to zero diameter
#     x = _sp.amin(_sp.r_[value[value!=-_sp.inf], -1e+9])
#     value[value==-_sp.inf] = x
#     value = value[pores]
#     return _sp.around(value, ndec)
#
#==============================================================================
