# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 13:33:11 2017

@author: islah
"""
from bwfpnm.Phases import GenericPhase


class Moisture(GenericPhase):
    r'''
    Creates Phase object with a default name 'moisture' and preset values for water
    P = 101325 Pa
    T = 293.15 K

    Parameters
    ----------
    network : bwfpnm Network object
        The network to which this phase object will be attached.

    Notes
    -----
    This explicit association is necessary so the Phase object can initialize
    data arrays of the correct size to store network data.
    The initial properties are all at std conditions of
    T = 293.15 K and P = 1 atm.

    Examples
    --------
    >>> import bwfpnm
    >>> pn = bwfpnm.Network.TestNet()
    >>> water = bwfpnm.Phases.Moisture(network=pn)
    '''
    def __init__(self, name='moisture', props=None, **kwargs):
        super().__init__(name=name, **kwargs)

        if props is None:
            props = self._default_params()

        for key in props:
            self['prop.'+key] = props[key]

    def _default_params(self):
        props = {'temperature': 20 + 273.15,            # K
                 'pressure': 101325.0,                  # Pa
                 'univ_gas_constant': 8.31447}          # J/mol.K
        return props


if __name__ =="__main__":
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    water = bwfpnm.Phases.Moisture(network=pn)
