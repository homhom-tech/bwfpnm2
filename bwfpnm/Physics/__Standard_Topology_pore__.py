# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:47:45 2015

@author: islah
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:24:52 2015

@author: islah
"""
from bwfpnm.Physics import models as pm
from bwfpnm.Physics import GenericPhysics

class Standard_Topology_pore(GenericPhysics):
    r"""
    Standard class for physics objects for Topological network model

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    pores and throats : array_like
        The pores and throats where this Physics object applies

    """

    def __init__(self, param='shapefactor', **kwargs):
        super().__init__(**kwargs)
        self._generate(param)

    def _generate(self, param):
        for phase in self._phases:
            temp = [item.split('.')[1] for item in phase.props()]
            if param in ['shapefactor', 'g', 'shape factor']:
                tmodel = pm.hydraulic_conductance.hp_shapefactor
                pc_method = 'g-based'
            elif param in ['radius', 'r', 'diameter']:
                tmodel = pm.hydraulic_conductance.hagen_poiseuille
                pc_method = 'r-based'
            if 'viscosity' in temp:
                self.models.add(propname='throat.hydraulic_conductance',
                                model=tmodel,
                                pore_area='throat.area')
                self.models.add(propname='throat.hydraulic_conductance_pore',
                                model=tmodel,
                                pore_area='pore.area')

            if 'surface_tension' in temp:
                self.models.add(propname='throat.capillary_pressure',
                                model=pm.capillary_pressure.washburn,
                                pore_diameter='throat.diameter',
                                method=pc_method)
                self.models.add(propname='pore.capillary_pressure',
                                model=pm.capillary_pressure.washburn,
                                pore_diameter='pore.diameter',
                                method=pc_method)


if __name__ == '__main__':
    print('none yet')
