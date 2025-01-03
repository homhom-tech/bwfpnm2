# -*- coding: utf-8 -*-
r"""
===============================================================================
Submodule -- vapor_pressure
===============================================================================

Methods for predicing the vapor pressure of pure species

"""
import scipy as sp

def pore(phase, Pc, with_RH=False, **kwargs):
    r"""
    Vapour pressure at pore as a function of Pc = f(rh)

    Parameters
    ----------


    """
    Pvsat = phase['prop.Pvsat']
    T = phase['prop.temperature']
    Rv = phase['prop.gas_constant']
    rho_l = phase['prop.water_density']

    RH = sp.exp(Pc/(rho_l*Rv*T))

    Pv = Pvsat*RH

    if with_RH:
        return Pv, RH
    else:
        return Pv

def saturation(phase, **kwargs):
    r"""Calculate Pvsat(T)
    """
    T = phase['prop.temperature'] - 273.15
    value = 288.68*(1.098 + T/100)**8.02
    return value

def saturation2(T, **kwargs):
    r"""Calculate Pvsat(T), T in Kelvin
    """
    T = T - 273.15
    value = 288.68*(1.098 + T/100)**8.02
    return value
