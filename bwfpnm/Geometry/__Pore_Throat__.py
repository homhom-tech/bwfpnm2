# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 15:00:51 2016

@author: islah

===============================================================================
Pore_Throat -- A modified 'stick & ball' geometrical model
===============================================================================

old version of PoreThroat: maintained for compatibility only

"""

from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry


class Pore_Throat(GenericGeometry):
    r"""
    Stick and Ball subclass of GenericGeometry.  This subclass is meant as a
    basic default geometry to get started quickly.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._generate(**kwargs)

    def _generate(self, psd_name='weibull_min', psd_shape=2.5, psd_loc=0,
                  psd_scale=0.5,
                  tsd_name='weibull_min', tsd_shape=2.5, tsd_loc=0,
                  tsd_scale=0.5, **kwargs):

        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        regen_mode='constant')
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        mode='min')

        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.sphere,
                        psd_name=psd_name,
                        psd_shape=psd_shape,
                        psd_loc=psd_loc,
                        psd_scale=psd_scale)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.sphere)

        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.cylinder,
                        tsd_name=tsd_name,
                        tsd_shape=tsd_shape,
                        tsd_loc=tsd_loc,
                        tsd_scale=tsd_scale)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.straight)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.cylinder)
        self.models.add(propname='throat.area',
                        model=gm.throat_area.cylinder)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.cylinder)
