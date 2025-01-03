# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 15:23:57 2015

@author: islah
"""
from OpenPNM.Geometry import GenericGeometry
import scipy as sp

class Test_MultiPoreThroat(GenericGeometry):
    r"""
    Default geometry for Pore-throat model with given geometry data.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, random=True, prandom=True,
                 tdiameter = 2e-6, tlength=2e-5,
                 macrorange=[-6, -3], microrange=[-9, -6], **kwargs):
        r"""
        Initialize
        """
        super().__init__(**kwargs)
        self._generate(random, prandom, tdiameter, tlength, macrorange,
                       microrange)

    def _generate(self, random, prandom, tdiameter, tlength,
                  macrorange, microrange, **kwargs):

        if prandom:
            pmacro = self._net['pore.macro']
            pmicro = self._net['pore.micro']
            Np = sp.sum(pmacro)
            lb, ub = macrorange
            pdiameter = 2*sp.logspace(lb, ub, Np)   # macro diameter
            sp.random.shuffle(pdiameter)
            pore_diameter = sp.zeros_like(self._net['pore.all'], dtype=float)
            pore_diameter[pmacro] = pdiameter

            mNp = sp.sum(pmicro)
            lb, ub = microrange
            mpdiameter = 2*sp.logspace(lb, ub, mNp)     # micro diameter
            sp.random.shuffle(mpdiameter)
            pore_diameter[pmicro] = mpdiameter
            self['pore.diameter'] = pore_diameter
        else:
            self['pore.diameter'] = sp.zeros_like(self._net['pore.all'], dtype=float)

        self['pore.area'] = sp.ones_like(self._net['pore.all'], dtype=float)
        self['pore.length'] = sp.ones_like(self._net['pore.all'], dtype=float)
        self['pore.volume'] = sp.zeros_like(self._net['pore.all'], dtype=float)

        if random:
            tmacro = self._net['throat.macro']
            tmicro = self._net['throat.micro']
            Nt = sp.sum(tmacro)
            tdiameter = 2*sp.logspace(-6, -3, Nt)
            sp.random.shuffle(tdiameter)
            tlength = 10*tdiameter

            mNt = sp.sum(tmicro)
            mtdiameter = 2*sp.logspace(-9, -6, mNt)
            sp.random.shuffle(mtdiameter)
            mtlength = 10*mtdiameter

        throat_diameter = sp.zeros_like(self._net['throat.all'], dtype=float)
        throat_length = sp.zeros_like(self._net['throat.all'], dtype=float)
        throat_diameter[tmacro] = tdiameter
        throat_diameter[tmicro] = mtdiameter
        throat_length[tmacro] = tlength
        throat_length[tmicro] = mtlength
#        self._net['throat.diameter'] = throat_diameter
        self['throat.diameter'] = throat_diameter
        self['throat.length'] = throat_length
        self['throat.area'] = sp.pi*self['throat.diameter']**2/4
        self['throat.volume'] = self['throat.length']*self['throat.area']




if __name__ == '__main__':
    #Run doc tests
#    import doctest
#    doctest.testmod(verbose=True)
    import bwfpnm
    pn = bwfpnm.Network.TestMultiNet(name='TestMultiNet', size=[2,2,2])
    Ps = pn.pores()
    Ts = pn.throats()
    geo = bwfpnm.Geometry.Test_MultiPoreThroat(network=pn, random=True, prandom=True,
                                           pores=Ps, throats=Ts,
                                           name='TestMultiThroat')
    bwfpnm.OpenPNM.Utilities.IO.VTK.save(pn, filename='TestMultiNetgeo2')
