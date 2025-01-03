# -*- coding: utf-8 -*-

from OpenPNM.Geometry import GenericGeometry
import scipy as sp

class Test_Throat(GenericGeometry):
    r"""
    Cylindrical throat geometry for throat model.
    Pore's diameter and volume are set to zero,
    while pore's area and length are set to one.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.
    random:     boolean.
                True, to enable random generation of throat's diameter and set its
                length to be 10x diameter. Otherwise, 'tdiameter' and 'tlength'
                have to be specified.

    """

    def __init__(self, random=False, trange=[-9, -3],
                 tdiameter = 2e-6, tlength=2e-5, **kwargs):
        r"""
        Initialize
        """
        super(Test_Throat, self).__init__(**kwargs)
        self._generate(random, tdiameter, tlength, trange[0], trange[1])

    def _generate(self, random, tdiameter, tlength, lb, ub, **kwargs):
        self['pore.diameter'] = 0
        self['pore.area'] = 1
        self['pore.length'] = 1
        self['pore.volume'] = 0
        if random:
            tdiameter = 2*sp.logspace(lb, ub, self.Nt)
            tlength = 10*tdiameter

        self['throat.diameter'] = tdiameter
        self['throat.length'] = tlength
        self['throat.area'] = sp.pi*self['throat.diameter']**2/4
        self['throat.volume'] = self['throat.length']*self['throat.area']


if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
