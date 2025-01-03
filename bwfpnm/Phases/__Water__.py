# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 13:48:24 2017

@author: islah
"""
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
    def __init__(self, name='water', props=None, **kwargs):
        super(Water, self).__init__(name=name, **kwargs)

        if props is None:
            props = self._default_params()

        props_none = []
        for key in props:
            if props[key] is not None:
                self['prop.'+key] = props[key]
#                self.props[key] = props[key]
            else:
                props_none.append(key)
        if props_none:
            self._generate(props_none)

    def _default_params(self):
        props = {'temperature': 20 + 273.15,            # K
                 'pressure': 101325.0,                  # Pa
                 'univ_gas_constant': 8.31447,          # J/mol.K
                 'molecular_weight': 0.01802,
                 'contact_angle': 0.0,
                 'critical_pressure': 2.2064e7,
                 'critical_temperature': 647.1,
                 'critical_volume': 0.003106,
                 'diffusivity': 1e-9}
        return props

    def _generate(self, props_none):
        if 'density' in props_none:
#            self.props['density'] = fm.density.water(self)
            self.models.add(propname='prop.density',
                            model=fm.density.water)              # kg/m3
        if 'molar_density' in props_none:
#            self.props['molar_density'] = fm.molar_density.standard(self)
            self.models.add(propname='prop.molar_density',
                            model=fm.molar_density.standard)     # mol/m3
        if 'surface_tension' in props_none:
#            self.props['surface_tension'] = fm.surface_tension.water(self)
            self.models.add(propname='prop.surface_tension',
                            model=fm.surface_tension.water)      # N/m
        if 'vapour_pressure' in props_none:
#            self.props['vapour_pressure'] = fm.vapor_pressure.antoine(self, A=10.1965,
#                            B=1730.63,
#                            C=-39.720)
            self.models.add(propname='prop.vapour_pressure',      # Pa
                            model=fm.vapor_pressure.antoine,
                            A=10.1965,
                            B=1730.63,
                            C=-39.720)
        if 'viscosity' in props_none:
#            self.props['viscosity'] = fm.viscosity.water(self)
            self.models.add(propname='prop.viscosity',
                            model=fm.viscosity.water)           # kg/m.s = Pa.s
##        self._molecular_weight = 0.01802             # kg/mol
##        self._contact_angle = 0.0                    # Degree
#        self['pore.molecular_weight'] = 0.01802             # kg/mol
#        self['pore.contact_angle'] = 0.0                    # Degree
#        self['pore.critical_pressure'] = 2.2064E7           # Pa
#        self['pore.critical_temperature'] = 647.1           # K
#        self['pore.critical_volume'] = 0.003106             # kg/m3
##        self._density = fm.density.water2(self._temperature)
##        self._molar_density = fm.molar_density.standard2(self._density,
##                                                         self._molecular_weight)
##        self._surface_tension = fm.surface_tension.water2(self._temperature)
##        self._viscosity = fm.viscosity.water2(self._temperature)
##        self._vapour_pressure = fm.vapor_pressure.antoine2(self._temperature,
##                                                          A=10.1965,
##                                                          B=1730.63,
##                                                          C=-39.720)
#        self['pore.diffusivity'] = 1e-9                     # m2/s

if __name__ == "__main__":
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    water = bwfpnm.Phases.Water(network=pn)
