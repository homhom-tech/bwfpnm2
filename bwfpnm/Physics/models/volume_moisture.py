# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:43:56 2015

@author: islah
"""
import scipy as _sp


def pvolume(phase,
            pc,
            pore_temperature='prop.temperature',
            pore_waterdensity='prop.water_density',
            pore_volume='pore.volume',
            pore_gasconstant='prop.gas_constant',
            pore_pvsat='prop.Pvsat',
            pore_occupancy='pore.occupancy_wp',
            **kwargs):
    r"""
    Calculate the volume of vapour in each pore at specified relative humidity.

    """
    V = phase._net[pore_volume]
    T = phase[pore_temperature]
    rho = phase[pore_waterdensity]
    Rv = phase[pore_gasconstant]
    Pvsat = phase[pore_pvsat]
    occupancy = phase[pore_occupancy]

    Va = V*occupancy
    rhoRT = rho*Rv*T
    rh = _sp.exp(pc/rhoRT)
    Pv = Pvsat*rh

    Mv = Va*Pv/(Rv*T)   # in kg
    Vv = Mv/rho         # in m^3
    return Vv


def tvolume(phase,
            pc,
            pore_temperature='prop.temperature',
            pore_waterdensity='prop.water_density',
            throat_volume='throat.volume',
            pore_gasconstant='prop.gas_constant',
            pore_pvsat='prop.Pvsat',
            throat_occupancy='throat.occupancy_wp',
            **kwargs):
    r"""
    Calculate the statistical thickness of each element at specified
    pc/relative humidity based on Bradley's equation.

    """
    V = phase._net[throat_volume]
    occupancy = phase[throat_occupancy]

    T = phase[pore_temperature]
    rho = phase[pore_waterdensity]
    Rv = phase[pore_gasconstant]
    Pvsat = phase[pore_pvsat]

    T = phase.interpolate_data(data=T)
    rho = phase.interpolate_data(data=rho)
    Rv = phase.interpolate_data(data=Rv)
    Pvsat = phase.interpolate_data(data=Pvsat)

    Va = V*occupancy
    rh = _sp.exp(pc/(rho*Rv*T))
    Pv = Pvsat*rh

    Mv = Va*Pv/(Rv*T)   # in kg
    Vv = Mv/rho         # in m^3
    return Vv
