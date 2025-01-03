# -*- coding: utf-8 -*-
from bwfpnm.Phases import GenericPhase
from bwfpnm.Phases import models as fm

class Air(GenericPhase):
    r"""
    Creates Phase object with preset models and values for air
    P = 101325 Pa
    T = 293.15 K
    contact angle = 0.0

    Parameters
    ----------
    network : bwfpnm Network object
        The network to which this phase object will be attached.

    Notes
    -----
    The initial properties are all at std conditions of T = 298 K and P = 1 atm.

    References
    ----------
    [1] E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity
    Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of Thermophysics,
    Vol. 25, No. 1, January 2004, pp. 21-69

    Examples
    --------
    >>> import bwfpnm
    >>> pn = bwfpnm.Network.TestNet()
    >>> air = bwfpnm.Phases.Air(network=pn)

    """
    def __init__(self,name=None,**kwargs):
        super(Air,self).__init__(name=name,**kwargs)

#        self._temperature = 20 + 273.15              # K
#        self._pressure = 101325.0                    # Pa
#        self._univ_gas_constant = 8.31447	          # J/mol.K
        self['prop.temperature'] = 20 + 273.15              # K
        self['prop.pressure'] = 101325.0                    # Pa
        self['prop.univ_gas_constant'] = 8.31447	    # J/mol.K
        self._generate()

    def _generate(self):
#        self._molecular_weight = 0.02896             # kg/mol
#        self._contact_angle = 0.0                    # Degree
#        self._gas_constant = 287			# J/kg.K
        self['prop.molecular_weight'] = 0.02896             # kg/mol
        self['prop.critical_pressure'] = 3.786E6            # Pa
        self['prop.critical_temperature'] = 132.5           # K
        self['prop.critical_volume'] = 0.002917             # kg/m3
        self['prop.contact_angle'] = 0.0                    # Degree
        self['prop.gas_constant'] = 287			# J/kg.K
#	self['prop.gas_constant'] = self['prop.univ_gas_constant']/self['prop.molecular_weight']	# kJ/kg.K = kPa.m^3/kg.K
        self.models.add(propname='prop.density',
                       model=fm.density.ideal_gas)          # kg/m3
        self.models.add(propname='prop.molar_density',
                       model=fm.molar_density.ideal_gas)    # mol/m3
        self['prop.diffusivity'] = 5.4785E-6                # m2/s
#        self.models.add(propname='prop.thermal_conductivity',# W/m.K
#                       model=fm.misc.polynomial,
#                       poreprop='prop.temperature',
#                       a=[0.00422791,0.0000789606,-1.56383E-08])
        self.models.add(propname='prop.viscosity',           # kg/m.s
                       model=fm.misc.polynomial,
                       poreprop='prop.temperature',
                       a=[0.00000182082,6.51815E-08,-3.48553E-11,1.11409E-14])

if __name__ =="__main__":
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    air = bwfpnm.Phases.Air(network=pn)
