# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 17:51:19 2015

@author: islah

===============================================================================
Submodule -- throat_diameter
===============================================================================

"""
import scipy.stats as spst
from scipy import power


def cylinder2(geometry, tsd_name, tsd_shape, tsd_loc, tsd_scale,
              throat_seed='throat.seed', tsd_offset=0, log_range=[-6, -3],
              **kwargs):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    """
    prob_fn = getattr(spst, tsd_name)
    P = prob_fn(tsd_shape, loc=tsd_loc, scale=tsd_scale)
    value = P.ppf(geometry[throat_seed]) + tsd_offset

    range_size = log_range[1]-log_range[0]
    value = value*range_size + log_range[0]
    value = power(10, value)/2
    return value

