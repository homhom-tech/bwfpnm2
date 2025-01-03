# -*- coding: utf-8 -*-

from OpenPNM.Geometry import GenericGeometry
import scipy as sp

class Test_PoreThroat(GenericGeometry):
    r"""
    Cylindrical geometry for general pore-throat model
    with random or specified pore's and throat's diameter and throat's length.

    The pore's length = pore's diameter.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.
    random:     boolean.
                True, to generate pore's radii spaced evenly on a log scale
                prange. Throat's diameters are half of the smallest radii of
                connected pores. Throat's length = 10x tdiameter.
                False: pdiameter, tdiameter, tlength must be supplied.
    prange:     Range of pore radii in log scale, [pmin, pmax].


    """

    def __init__(self, random=False, prange=[-9, -3],
                 pdiameter=2e-5, tdiameter=2e-6, tlength=1e-4, **kwargs):
        r"""
        Initialize
        """
        super().__init__(**kwargs)
        self._generate(random, pdiameter, tdiameter, tlength,
                       prange[0], prange[1])


    def _generate(self, random, pdiameter, tdiameter, tlength,
                  ul, ub, **kwargs):
        if random:
            # generate different pore diameter between 2e-9 and 2e-3
            pdiameter = 2*sp.logspace(ul, ub, self.Np)
            # throat diameter = min(diameter(connected pores))/2
            tind = self._net.throats()
            conns = self._net.find_connected_pores(tind)
            tdiameter = sp.amin(pdiameter[conns], axis=1)/2
            tlength = 10*tdiameter

        self['pore.diameter'] = pdiameter
        self['pore.area'] = sp.pi*pdiameter**2/4
        self['pore.length'] = pdiameter
        self['pore.volume'] = self['pore.length']*self['pore.area']

        self['throat.diameter'] = tdiameter
        self['throat.length'] = tlength

        self['throat.area'] = sp.pi*tdiameter**2/4
        self['throat.volume'] = tlength*self['throat.area']


if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
