"""
module Physics
===============================================================================

"""

from bwfpnm.Physics import models as pm
from bwfpnm.Physics import GenericPhysics

class Standard(GenericPhysics):
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

    def __init__(self,**kwargs):
        super(Standard,self).__init__(**kwargs)
        self._generate()

    def _generate(self):
        for phase in self._phases:
            temp = [item.split('.')[1] for item in phase.props()]
            if 'viscosity' in temp:
                self.models.add(propname='throat.hydraulic_conductance',
                                model=pm.hydraulic_conductance.thagen_poiseuille)
                self.models.add(propname='pore.hydraulic_conductance',
                                model=pm.hydraulic_conductance.phagen_poiseuille)
#            if 'diffusivity' in temp:
#                self.models.add(propname='throat.diffusive_conductance',
#                                model=pm.diffusive_conductance.tbulk_diffusion)
#                self.models.add(propname='pore.diffusive_conductance',
#                                model=pm.diffusive_conductance.pbulk_diffusion)
            if 'surface_tension' in temp:
                self.models.add(propname='throat.capillary_pressure',
                                model=pm.capillary_pressure.twashburn)
                self.models.add(propname='pore.capillary_pressure',
                                model=pm.capillary_pressure.pwashburn)


if __name__ == '__main__':
    print('none yet')
