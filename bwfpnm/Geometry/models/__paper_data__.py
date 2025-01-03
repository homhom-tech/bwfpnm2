# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 09:42:37 2014

@author: islah
"""


class material(object):
    pass


def ceramicbrick_carmeliet01(case='Wetting'):
    """Data based on Carmeliet & Roels 2001:
    'Determination of the isothermal moisture transport properties
    of porous building materials'

    argument:   case = ['Wetting', 'Drainage']
    return:     Ceramicbrick(an object)
    """
    Ceramicbrick = material()
    Ceramicbrick.name = 'Ceramicbrick'
    Ceramicbrick.case = case
    # --- basic moisture properties ---
    Ceramicbrick.open_porosity = 0.24
    Ceramicbrick.bulk_density = 2005
    Ceramicbrick.w_sat = 240
    Ceramicbrick.w_cap = 160
    Ceramicbrick.A_cap = 0.165
    Ceramicbrick.K_sat = 1.3e-7
    if case == 'Wetting':
        # --- network parameteres: Wetting ---
        Ceramicbrick.vapor_resistance = [48.3, 31.2, 8.7]
        Ceramicbrick.radius_max = 1.58e-5
        Ceramicbrick.radius_min = 3.16e-9
        Ceramicbrick.radius_r1 = 7.08e-8
        Ceramicbrick.number_of_scales = 5
        Ceramicbrick.pore_shape_ratio = -1
        # --- capillary presssure curve: Wetting ---
        Ceramicbrick.number_of_distribution = 3
        Ceramicbrick.parameter_c = [2.763e-5, 1.715e-5, 3.835e-4]
        Ceramicbrick.exponent_n = [1.691, 4.457, 1.320]
        Ceramicbrick.weight_factor = [0.4544, 0.4913, 0.0543]
    elif case == 'Drainage':
        # --- network parameteres: Drainage ---
        Ceramicbrick.vapor_resistance = [0.29, 0.68, 0.9]
        Ceramicbrick.radius_max = 1.58e-5
        Ceramicbrick.radius_min = 3.16e-9
        Ceramicbrick.radius_r1 = 6.31e-6
        Ceramicbrick.number_of_scales = 11
        Ceramicbrick.pore_shape_ratio = -1
        # --- capillary presssure curve: Drainage ---
        Ceramicbrick.number_of_distribution = 3
        Ceramicbrick.parameter_c = [9.9e-6, 2.4e-5, 3.84e-4]
        Ceramicbrick.exponent_n = [1.638, 4.0, 1.5]
        Ceramicbrick.weight_factor = [0.165, 0.7808, 0.0542]
    return Ceramicbrick


def calciumsilicate_carmeliet01(case='Wetting'):
    """Data based on Carmeliet & Roels 2001:
    'Determination of the isothermal moisture transport properties
    of porous building materials'

    argument: case = ['Wetting', 'Drainage']
    return: Calciumsilicate(an object)
    """
    Calciumsilicate = material()
    Calciumsilicate.name = 'Calciumsilicate'
    Calciumsilicate.case = case
    # --- basic moisture properties ---
    Calciumsilicate.open_porosity = 0.319
    Calciumsilicate.bulk_density = 1803
    Calciumsilicate.w_sat = 319
    Calciumsilicate.w_cap = 239
    Calciumsilicate.A_cap = 0.043
    Calciumsilicate.K_sat = 7.9e-8
    if case == 'Wetting':
        # --- network parameteres: Wetting ---
        Calciumsilicate.vapor_resistance = [23.0, 12.0, 3.0]
        Calciumsilicate.radius_max = 1.0e-4
        Calciumsilicate.radius_min = 4.47e-10
        Calciumsilicate.radius_r1 = 2.51e-8
        Calciumsilicate.number_of_scales = 7
        Calciumsilicate.pore_shape_ratio = -1.4
        # --- capillary presssure curve: Wetting ---
        Calciumsilicate.number_of_distribution = 3
        Calciumsilicate.parameter_c = [3.4e-7, 2.5e-6, 4.5e-5]
        Calciumsilicate.exponent_n = [1.85, 1.6, 4.0]
        Calciumsilicate.weight_factor = [0.274, 0.653, 0.073]
    elif case == 'Drainage':
        # --- network parameteres: Drainage ---
        Calciumsilicate.vapor_resistance = [0.3, 0.7, 0.92]
        Calciumsilicate.radius_max = 1.0e-4
        Calciumsilicate.radius_min = 4.47e-10
        Calciumsilicate.radius_r1 = 2.51e-7
        Calciumsilicate.number_of_scales = 11
        Calciumsilicate.pore_shape_ratio = -1.4
        # --- capillary presssure curve: Drainage ---
        Calciumsilicate.number_of_distribution = 3
        Calciumsilicate.parameter_c = [3.4e-7, 2.5e-6, 8.0e-5]
        Calciumsilicate.exponent_n = [1.85, 1.6, 4.0]
        Calciumsilicate.weight_factor = [0.21, 0.5, 0.29]
    return Calciumsilicate


def calciumsilicate_carmeliet99():
    """Data based on Carmeliet 1999:
    'A multiscale network model for simulating moisture transfer properties of
     porous building materials'

    argument: -
    return: Calciumsilicate(an object)
    """
    Calciumsilicate = material()
    Calciumsilicate.name = 'Calciumsilicate'
    # --- Distribution Parameters ---
    Calciumsilicate.number_of_distribution = 3
    Calciumsilicate.parameter_c = [3.55e-7, 1.19e-5, 1.15e-4]
    Calciumsilicate.exponent_n = [2.026, 1.238, 6.299]
    Calciumsilicate.weight_factor = [0.27, 0.46, 0.27]
    Calciumsilicate.w_sat = 359

    # --- Multiscale network model parameters ---
    Calciumsilicate.net_size_s = 30
    Calciumsilicate.sim_num_n = 40
    Calciumsilicate.pore_shape_ratio = -1
    Calciumsilicate.connct_z = 8
    Calciumsilicate.radius_max = 48.0e-6
    Calciumsilicate.radius_min = 1.0e-9

    # --- Basic moisture properties ---
    Calciumsilicate.bulk_density = 1689
    Calciumsilicate.open_porosity = 0.359
    Calciumsilicate.K_sat = 7.9e-8
    Calciumsilicate.w_cap = 260

    return Calciumsilicate


def ceramicbrick_carmeliet99():
    """Data based on Carmeliet 1999:
    'A multiscale network model for simulating moisture transfer properties of
     porous building materials'

    argument: -
    return: Calciumsilicate(an object)
    """
    Ceramicbrick = material()
    Ceramicbrick.name = 'Ceramicbrick'
    # --- Distribution Parameters ---
    Ceramicbrick.number_of_distribution = 2
    Ceramicbrick.parameter_c = [9.77e-6, 2.0e-5]
    Ceramicbrick.exponent_n = [1.78, 5.5]
    Ceramicbrick.weight_factor = [0.25, 0.75]      # bi-modal
#    Ceramicbrick.weight_factor = [0., 1.]       # uni-modal --> single scale
    Ceramicbrick.w_sat = 242

    # --- Multiscale network model parameters ---
    Ceramicbrick.net_size_s = 30
    Ceramicbrick.sim_num_n = 40
    Ceramicbrick.pore_shape_ratio = -1
    Ceramicbrick.connct_z = 8
    Ceramicbrick.radius_max = 15.8e-6
    Ceramicbrick.radius_min = 1.0e-9

    # --- Basic moisture properties ---
    Ceramicbrick.bulk_density = 2005
    Ceramicbrick.open_porosity = 0.24
    Ceramicbrick.K_sat = 1.3e-7
    Ceramicbrick.w_cap = 157

    return Ceramicbrick
