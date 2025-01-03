r"""
===============================================================================
Submodule -- diffusivity
===============================================================================

"""
#import scipy as sp

def vapour_air(phase, **kwargs):
    r"""
    Vapour permeability of air

    Hans:
    \delta_{v,a} = 2.262/(R_v*T*P)*(T/273.15)**1.81

    For T = 20C, value = 1.9e-10 kg/msPa

    Here, the denominators (R_v*T) are excluded. They will be appeared in
    the main equation.

    Parameters
    ----------

    """
    T = phase['prop.temperature']
    P = phase['prop.pressure']

    value = 2.262/P*(T/273.15)**1.81            # m2/s <- checked!
    return value

def vapour_air2(T, P, **kwargs):
    r"""
    Vapour permeability of air

    Hans:
    \delta_{v,a} = 2.262/(R_v*T*P)*(T/273.15)**1.81

    For T = 20C, value = 1.9e-10 kg/msPa

    Here, the denominators (R_v*T) are excluded. They will be appeared in
    the main equation.

    Parameters
    ----------

    """
    value = 2.262/P*(T/273.15)**1.81            # m2/s <- checked!
    return value