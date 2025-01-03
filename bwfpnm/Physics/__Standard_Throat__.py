# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 00:19:16 2015

@author: islah
"""

from bwfpnm.Physics import models as pm
from bwfpnm.Physics import GenericPhysics

class Standard_Throat(GenericPhysics):
    r"""
    Base class to generate a generic Physics object.  The user must specify models
    and parameters for the all the properties they require. Classes for several
    common Physics are included with OpenPNM and can be found under OpenPNM.Physics.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    pores and throats : array_like
        The pores and throats where this Physics object applies

    """

    def __init__(self, **kwargs):
        super(Standard_Throat, self).__init__(**kwargs)
        self._generate()

    def _generate(self):
        for phase in self._phases:
            temp = [item.split('.')[1] for item in phase.props()]
            if 'viscosity' in temp:
                self.models.add(propname='throat.hydraulic_conductance',
                                model=pm.hydraulic_conductance.thagen_poiseuille)

            if 'surface_tension' in temp:
                self.models.add(propname='throat.capillary_pressure',
                                model=pm.capillary_pressure.twashburn)


if __name__ == '__main__':
    print('none yet')
