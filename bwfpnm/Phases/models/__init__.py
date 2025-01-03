r"""
*******************************************************************************
models -- Functions for calculating thermofluid properties
*******************************************************************************

Contents
--------
This module contains methods for estimating phase properties

"""
from OpenPNM.Phases.models import (density, molar_density, misc,
                                   vapor_pressure, viscosity)
from . import diffusivity
from . import surface_tension
from . import vapour_pressure

