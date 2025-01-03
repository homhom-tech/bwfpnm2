# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:16:02 2015

@author: islah

"""
from bwfpnm.Physics import models as pm
from bwfpnm.Physics import GenericPhysics

class Standard_Topology_pore_eq(GenericPhysics):
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
            if 'viscosity' in temp:
                if param=='shapefactor':
                    tmodel = pm.hydraulic_conductance.thp_shapefactor
                    tmodelp = pm.hydraulic_conductance.thp_shapefactor_pore
                elif param == 'radius':
                    tmodel = pm.hydraulic_conductance.thagen_poiseuille

                self.models.add(propname='throat.hydraulic_conductance',
                                model=tmodel,
                                throat_area='throat.area_eq')
                self.models.add(propname='throat.hydraulic_conductance_pore',
                                model=tmodelp,
                                pore_area='pore.area_eq')

            if 'surface_tension' in temp:
                self.models.add(propname='throat.capillary_pressure',
                                model=pm.capillary_pressure.twashburn_eq)
                self.models.add(propname='pore.capillary_pressure',
                                model=pm.capillary_pressure.pwashburn_eq)


if __name__ == '__main__':
    print('none yet')
