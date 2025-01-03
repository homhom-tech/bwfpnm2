# -*- coding: utf-8 -*-
from . import GenericPhase
from . import models as fm


class Water(GenericPhase):
    r'''
    Creates Phase object with a default name 'water' and preset values for water
    P = 101325 Pa
    T = 293.15 K
    contact angle = 0.0

    Parameters
    ----------
    network : bwfpnm Network object
        The network to which this phase object will be attached.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    Examples
    --------
    >>> import bwfpnm
    >>> pn = bwfpnm.Network.TestNet()
    >>> water = bwfpnm.Phases.Water(network=pn)
    '''
    def __init__(self, name='water', **kwargs):
        super(Water, self).__init__(name=name, **kwargs)

#        self._temperature = 20 + 273.15              # K
#        self._pressure = 101325.0                    # Pa
        self['pore.temperature'] = 20 + 273.15              # K
        self['pore.pressure'] = 101325.0                    # Pa
        self._generate()

    def _generate(self):
#        self._molecular_weight = 0.01802             # kg/mol
#        self._contact_angle = 0.0                    # Degree
        self['pore.molecular_weight'] = 0.01802             # kg/mol
        self['pore.contact_angle'] = 0.0                    # Degree
        self['pore.critical_pressure'] = 2.2064E7           # Pa
        self['pore.critical_temperature'] = 647.1           # K
        self['pore.critical_volume'] = 0.003106             # kg/m3
#        self._density = fm.density.water2(self._temperature)
#        self._molar_density = fm.molar_density.standard2(self._density,
#                                                         self._molecular_weight)
#        self._surface_tension = fm.surface_tension.water2(self._temperature)
#        self._viscosity = fm.viscosity.water2(self._temperature)
#        self._vapour_pressure = fm.vapor_pressure.antoine2(self._temperature,
#                                                          A=10.1965,
#                                                          B=1730.63,
#                                                          C=-39.720)
        self.models.add(propname='pore.density',
                        model=fm.density.water)              # kg/m3
        self.models.add(propname='pore.molar_density',
                        model=fm.molar_density.standard)     # mol/m3
        self.models.add(propname='pore.surface_tension',
                        model=fm.surface_tension.water)      # N/m
        self.models.add(propname='pore.vapour_pressure',      # Pa
                        model=fm.vapor_pressure.antoine,
                        A=10.1965,
                        B=1730.63,
                        C=-39.720)
        self.models.add(propname='pore.viscosity',
                        model=fm.viscosity.water)            # kg/m.s = Pa.s
        self['pore.diffusivity'] = 1e-9                     # m2/s

if __name__ == "__main__":
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    water = bwfpnm.Phases.Water(network=pn)
