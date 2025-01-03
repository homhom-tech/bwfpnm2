# -*- coding: utf-8 -*-
from bwfpnm.Phases import GenericPhase
from bwfpnm.Phases import models as fm
from bwfpnm.Phases import models as fmb

class Vapour(GenericPhase):
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
    def __init__(self, name='vapour', props=None, **kwargs):
        super(Vapour,self).__init__(name=name, **kwargs)

        if props is None:
            props = self._default_params()

        props_none = []
        for key in props:
            if props[key] is not None:
#                self['pore.'+key] = props[key]
                self['prop.'+key] = props[key]
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
                 'diffusivity': 2.54e-5, # m2/s 1e-9,
                 'density': 0.0173,
                 'viscosity': 1.3e-5,
                 'Pvsat': 2338,
                 'gas_constant': 461.4}
        return props

    def _generate(self, props_none):
        r"""
        Refs:
        - Y. Cengel: Appendix
        """
        if 'Pvsat' in props_none:
            self.models.add(propname='prop.Pvsat',
                            model=fmb.vapour_pressure.saturation)   # P
        if 'molar_density' in props_none:
            self.models.add(propname='prop.molar_density',
                            model=fm.molar_density.standard)     # mol/m3
        if 'diffusivity' in props_none:
            self.models.add(propname='prop.diffusivity',
                            model=fmb.diffusivity.vapour_air)       # m2/s
        if 'permeability' in props_none:
            DAB = self['prop.diffusivity']
            R_v = self['prop.gas_constant']
            T = self['prop.temperature']
            self['prop.permeability'] = DAB/(R_v*T)            # s <- checked

#        self._molecular_weight = 0.01802             # kg/mol
#        self['pore.molecular_weight'] = 0.01802             # kg/mol
#        Ru = self['pore.univ_gas_constant']
#        M = self['pore.molecular_weight']
#        Ru = self._univ_gas_constant
#        M = self._molecular_weight
#        self._gas_constant = Ru/M	          # kJ/kg.K = kPa.m^3/kg.K
#        self._density = 0.0173			     # kg/m^3
#        self._contact_angle = 0.0               # Degree
#        self['pore.gas_constant'] = Ru/M	          # kJ/kg.K = kPa.m^3/kg.K
#        self['pore.density'] = 0.0173			     # kg/m^3
#        self['pore.contact_angle'] = 0.0                    # Degree

#        self._Pvsat = fmb.vapour_pressure.saturation2(self._temperature)
#        self._molar_density = fm.molar_density.standard2(self._density,
#                                                         self._molecular_weight)
#        self._diffusivity = fmb.diffusivity.vapour_air2(self._temperature,
#                                                        self._pressure)
#        DAB = self._diffusivity
#        R_v = self._gas_constant
#        T = self._temperature
#        self._permeability = DAB/(R_v*T)            # s <- checked
#        self._viscosity = 1.3e-5		          # kg/m.s


#        self['pore.viscosity'] = 1.3e-5		          # kg/m.s
#        self['pore.Pvsat'] = 2339			          # P
#        self['pore.diffusivity'] = 1e-9
#        self['pore.gas_constant'] = 461.4		# J/kg.K



if __name__ =="__main__":
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    water = bwfpnm.Phases.Water(network=pn)
