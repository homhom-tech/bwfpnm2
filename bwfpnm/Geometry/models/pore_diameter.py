# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 17:18:14 2015

@author: islah

===============================================================================
pore_diameter
===============================================================================

"""
import scipy.stats as spst
from scipy import power

def sphere2(geometry, psd_name, psd_shape, psd_loc, psd_scale,
           pore_seed='pore.seed', psd_offset=0, log_range=[-6, -3],
           **kwargs):
    r"""
    Calculate pore diameter from given seed values.

    Parameters
    ----------
    geometry : OpenPNM Geometry Object
        The Geometry object which this model is associated with. This controls
        the length of the calculated array, and also provides access to other
        necessary geometric properties.

    psd_name : string
        The name of the statistical distribution to use. This model uses the
        Scipy.stats module, so any of the distributions available there are
        suitable options.

    psd_shape, loc and scale : scalars
        The parameters to send to the selected statistics model.  Most of the
        Scipy.stats models accept the same keyword arguments.  Note that the
        psd_ prefix is added by OpenPNM to indicate 'pore size distribution'.

    psd_offset : scalar
        Controls the minimum value in the pore size distribution by shifting
        the entire set of values by the given offset.  This is useful for
        avoiding pore sizes too close to zero.

    log_range : list/array of 2 components
        The min and max log values

    """
    prob_fn = getattr(spst, psd_name)
    P = prob_fn(psd_shape, loc=psd_loc, scale=psd_scale)
    value = P.ppf(geometry[pore_seed]) + psd_offset

    range_size = log_range[1]-log_range[0]
    value = value*range_size + log_range[0]
    value = power(10, value)
    return value
