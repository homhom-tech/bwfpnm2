# -*- coding: utf-8 -*-
from OpenPNM.Phases import GenericPhase
from OpenPNM.Phases import models as fm


class Mercury(GenericPhase):
    r"""
    Creates Phase object with a default name 'Hg' and preset values for
    mercury.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    References
    ----------
    [1] Thermophysical Properties of Materials for Nuclear Engineering: IAEA,
        Vienna, 2008. ISBN 978-92-0-106508-7:

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> hg = OpenPNM.Phases.Mercury(network=pn)

    """
    def __init__(self, name=None, props=None, **kwargs):
        super().__init__(name=name, **kwargs)
        if props is None:
            props = self._default_params()
        props_none = []
        for key in props:
            if props[key] is not None:
                self['prop.'+key] = props[key]
#                self.props[key] = props[key]
            else:
                props_none.append(key)
#        if props_none:
#            self._generate(props_none)
        self._generate()

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

    def _generate(self):
        self['pore.molecular_weight'] = 0.2006                         # kg/mol
        self['pore.critical_pressure'] = 1.662E8                       # Pa
        self['pore.critical_temperature'] = 1733                       # K
        self['pore.critical_volume'] = 0.000189                        # kg/m3
        self['pore.contact_angle'] = 140                               # Degree
        self.models.add(propname='prop.vapor_pressure',                 # Pa
                        model=fm.vapor_pressure.antoine,
                        A=9.85767, B=3007.129, C=-10.001)
        self.models.add(propname='prop.density',                        # kg/m3
                        model=fm.misc.linear,
                        poreprop='prop.temperature',
                        b=14280.9, m=-2.47004)
        self.models.add(propname='prop.molar_density',
                        model=fm.molar_density.standard)                # mol/m3
        self.models.add(propname='prop.surface_tension',                # N/m
                        model=fm.misc.linear,
                        poreprop='prop.temperature',
                        b=0.56254, m=-0.00028)
        self.models.add(propname='prop.thermal_conductivity',           # W/m.K
                        model=fm.misc.polynomial,
                        poreprop='prop.temperature',
                        a=[3.98691, 0.0170967, -0.0000063862])
        self.models.add(propname='prop.viscosity',                      # kg/m.s
                        model=fm.misc.polynomial,
                        poreprop='prop.temperature',
                        a=[0.00355837, -0.0000100131, 1.23684E-08, -5.16836E-12])
