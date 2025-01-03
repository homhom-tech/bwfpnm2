r"""
===============================================================================
Submodule -- surface_tension
===============================================================================

Some text here?

"""

import scipy as sp

def water(phase,**kwargs):
    r"""
    Calculates surface tension = f(T)
    """
    T = phase['prop.temperature']
    value = (122.3 - 0.17*T)*1e-3	# N/m
    return value

def water2(T,**kwargs):
    r"""
    Calculates surface tension = f(T)
    """
    value = (122.3 - 0.17*T)*1e-3	# N/m
    return value




