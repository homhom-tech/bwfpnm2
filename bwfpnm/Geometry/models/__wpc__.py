# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 09:20:59 2014

@author: islah

Method: Moisture Retention Curve
"""
#import os
#import sys
#sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import numpy as np

__author__ = '\n'.join(['Muhammad Islahuddin '
                        '<islah.islahuddin@bwk.kuleuven.be>'])

__all__ = ['pc_from_r', 'r_from_pc', 'darcy_q', 'k_liq', 'k_vap',
           'rh_from_pc', 'pv_sat', 'penet_depth']


def pc_from_r(rad, surftension=7.25e-2, angle=0):
    """ Capillary Pressure Pc [Pa] --> Young-Laplace eq with Pc = Pl - Pg """
    return -2*surftension*np.cos(angle)/np.array(rad)


def r_from_pc(pc, surftension=7.25e-2, angle=0):
    """Young-Laplace eq"""
    return -2*surftension*np.cos(angle)/np.array(pc)


def darcy_q(IntPermeability, dP, dx, DynViscWtr=1.002e-3):
    """ q [m/s] = -(k/miu)*dP/dx, k = intrinsic permeability """
    return -IntPermeability/DynViscWtr*dP/dx


def k_liq(RadPore, DensityWater=1e+3, DynViscWtr=1.002e-3):
    """Liquid permeability [s] from Hagen-Poiseuille eq """
    return DensityWater*RadPore**2/(8*DynViscWtr)


def k_vap(Pv, VapPermeabiltyInAir=1.9e-10, Temp=293.15, DensityWater=1e+3,
          Rv=462):
    """ Vapour permeability [s] = da*rho_v/rho_l """
    return VapPermeabiltyInAir*Pv/(DensityWater*Rv*Temp)


def rh_from_pc(pc, Temp=293.15, DensityWater=1000, Rv=462):
    """ Relative humidity [-] from Kelvin equation, RH_star for a pore
        T in K = 20+273.15 K """
    return np.exp(np.array(pc)/DensityWater/Rv/Temp)


def pv_sat(Temp=20):
    """ y= Pvsat [Pa], Temp in C """
    return 288.68*(1.098 + Temp/100)**8.02


def penet_depth(rad, DynVisc=1.002e-3, surftension=7.25e-2, angle=0):
    """ Water penetration depth B [m/s^0.5] """
    return np.sqrt(rad*surftension*np.cos(angle)/DynVisc/2)


if __name__ == "__main__":
    print('This is the wpc module to be imported from main code')
