# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:29:02 2016

@author: islah
"""

from OpenPNM.Geometry import models as gm
from bwfpnm.Geometry import models as gmb
from OpenPNM.Geometry import GenericGeometry
from bwfpnm.Geometry import PoreThroat
import scipy.stats as spst
from scipy import power
import scipy as _sp


class MultiPoreThroat(GenericGeometry):
    r"""
    Default geometry for Pore-throat model with given geometry data.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    """

    def __init__(self, **kwargs):
        r"""
        Initialize
        """
        super().__init__(**kwargs)
        self._generate(**kwargs)

    def _generate(self, macrorange=[-6, -4], microrange=[-9, -7], **kwargs):
        pn = self._net
        pmacro, pmicro = pn.pores('macro'), pn.pores('micro')
        tmacro, tmicro = pn.throats('macro'), pn.throats('micro')

        props = ['pore.seed', 'pore.diameter', 'pore.area', 'pore.volume',
                 'throat.seed', 'throat.diameter',
                 'throat.volume', 'throat.area']
        macro = {'pore': pmacro, 'throat': tmacro}
        micro = {'pore': pmicro, 'throat': tmicro}

        macro_prop = self._subscale(pmacro, tmacro, log_range=macrorange)
        micro_prop = self._subscale(pmicro, tmicro, log_range=microrange)

        self._create_props(props, macro, micro, macro_prop, micro_prop)

        tlength = self._tlength(self._net, pn.throats(),
                                pore_diameter='pore.diameter', L_negative=1e-9)
        tsurfarea = _sp.constants.pi*self['throat.diameter']*tlength

        self['throat.length'] = tlength[0]
        self['throat.surface_area'] = tsurfarea

        # p/t.shapefactor, p/t.shapefactor_constant,
        # throat.porelengths, throat.lengthtotal
        self.add_model(propname='pore.shapefactor',
                       model=gmb.shapefactor.shapefactor,
                       radius=self['pore.diameter']/2)
        self.add_model(propname='throat.shapefactor',
                       model=gmb.shapefactor.shapefactor,
                       radius=self['throat.diameter']/2)
        self.add_model(propname='pore.shapefactor_constant',
                       model=gmb.shapefactor.shapefactor,
                       radius=self['pore.shapefactor'])
        self.add_model(propname='throat.shapefactor_constant',
                       model=gmb.shapefactor.shapefactor,
                       radius=self['throat.shapefactor'])
        self['throat.porelengths'] = self._tporelength(self._net)
        self['throat.lengthtotal'] = tlength[1]

    def _create_props(self, props, macro, micro, macro_prop, micro_prop):
        for prop in props:
            self[prop] = 0.
            self[prop] = 0.
            split = prop.split('.')
            self[prop][macro[split[0]]] = macro_prop[prop]
            self[prop][micro[split[0]]] = micro_prop[prop]

    def _subscale(self, pores, throats, log_range=[-6, -3],
                  psd_name='weibull_min', psd_shape=2.5, psd_loc=0,
                  psd_scale=0.5, psd_offset=0,
                  tsd_name='weibull_min', tsd_shape=2.5, tsd_loc=0,
                  tsd_scale=0.5, tsd_offset=0, **kwargs):

        np, nt = _sp.size(pores), _sp.size(throats)
        pseed = self._random(np)
        tseed = self._neighbor(throats, pores, pseed, mode='min')
        pdia = self._diameter(np, psd_name, psd_shape, psd_loc, psd_scale,
                         pseed, psd_offset, log_range)
        tdia = self._diameter(nt, tsd_name, tsd_shape, tsd_loc, tsd_scale,
                         tseed, tsd_offset, log_range)
        parea = _sp.constants.pi/4*(pdia)**2
        tarea = _sp.constants.pi/4*(tdia)**2
        pvol = _sp.pi/6*pdia**3
        tvol = _sp.pi/6*tdia**3

        props = {'pore.seed': pseed, 'pore.diameter': pdia, 'pore.area': parea,
                 'pore.volume': pvol,
                 'throat.seed': tseed, 'throat.diameter': tdia,
                 'throat.area': tarea, 'throat.volume': tvol}
        return props

    def _diameter(self, ndia, psd_name, psd_shape, psd_loc, psd_scale,
                  pore_seed, psd_offset=0, log_range=[-6, -3],
                  **kwargs):
        r"""
        Calculate pore diameter from given seed values.

        Parameters
        ----------
        geometry : OpenPNM Geometry Object
            The Geometry object which this model is associated with. This controls
            the length of the calculated array, and also provides access to other
            necessary geometric properties.

        psd_name : string
            The name of the statistical distribution to use. This model uses the
            Scipy.stats module, so any of the distributions available there are
            suitable options.

        psd_shape, loc and scale : scalars
            The parameters to send to the selected statistics model.  Most of the
            Scipy.stats models accept the same keyword arguments.  Note that the
            psd_ prefix is added by OpenPNM to indicate 'pore size distribution'.

        psd_offset : scalar
            Controls the minimum value in the pore size distribution by shifting
            the entire set of values by the given offset.  This is useful for
            avoiding pore sizes too close to zero.

        log_range : list/array of 2 components
            The min and max log values

        """
        prob_fn = getattr(spst, psd_name)
        P = prob_fn(psd_shape, loc=psd_loc, scale=psd_scale)
        value = P.ppf(pore_seed) + psd_offset

        range_size = log_range[1]-log_range[0]
        value = value*range_size + log_range[0]
        value = power(10, value)
        return value

    def _random(self, number, seed=None, num_range=[0, 1], **kwargs):
        r"""
        Assign random number to pore bodies
        note: should this be called 'poisson'?
        """
        range_size = num_range[1]-num_range[0]
        range_min = num_range[0]
        _sp.random.seed(seed=seed)
        value = _sp.random.rand(number,)
        value = value*range_size + range_min
        return value

    def _neighbor(self, throats, pores, pseed, mode='min', **kwargs):
        r"""
        Adopt a value based on the neighboring pores
        """
        P12 = self._net.find_connected_pores(throats)
        temp = _sp.zeros(pores.max()+1)
        temp[pores] = pseed
        pvalues = temp[P12]

#        pvalues = self._net[pore_prop][P12]
        if mode == 'min':
            value = _sp.amin(pvalues, axis=1)
        elif mode == 'max':
            value = _sp.amax(pvalues, axis=1)
        elif mode == 'mean':
            value = _sp.mean(pvalues, axis=1)
        return value

    def _tlength(self, network, throats, pore_diameter='pore.diameter',
                 L_negative=1e-9, **kwargs):
        r"""
        Calculate throat length

        Parameters
        ----------
        L_negative : float
            The default throat length to use when negative lengths are found.  The
            default is 1 nm.  To accept negative throat lengths, set this value to
            ``None``.
        """
        # Initialize throat_property['length']
        pore1 = network['throat.conns'][:, 0]
        pore2 = network['throat.conns'][:, 1]
        C1 = network['pore.coords'][pore1]
        C2 = network['pore.coords'][pore2]
        E = _sp.sqrt(_sp.sum((C1-C2)**2, axis=1))  # Euclidean distance between pores
        D1 = network[pore_diameter][pore1]
        D2 = network[pore_diameter][pore2]
        value = E-(D1+D2)/2.
        value = value[throats]
        if _sp.any(value < 0) and L_negative is not None:
            print('Negative throat lengths are calculated. Arbitrary positive \
                   length assigned: ' + str(L_negative))
            Ts = _sp.where(value < 0)[0]
            value[Ts] = L_negative
        return (value, E)

    def _tporelength(self, network, pore_diameter='pore.diameter',
                     L_negative=1e-6, **kwargs):
        r"""
        Calculate throat's pore lengths (radii of each pore ends)

        Parameters
        ----------
        L_negative : float
            The default throat length to use when negative lengths are found.  The
            default is 1 micron.  To accept negative throat lengths, set this value to
            ``None``.
        """
        # Initialize throat_property['length']
        conns = network['throat.conns']
        value = network[pore_diameter][conns]/2

        if _sp.any(value < 0) and L_negative is not None:
            print('Negative throat lengths are calculated. Arbitrary positive \
                   length assigned: ' + str(L_negative))
            Ts = _sp.where(value < 0)[0]
            value[Ts] = L_negative
        return value

if __name__ == '__main__':
    #Run doc tests
#    import doctest
#    doctest.testmod(verbose=True)
    import bwfpnm
    pn = bwfpnm.Network.MultiNet(name='MultiNet', size=[2,2,2], size_micro=[3,3,3])
    Ps = pn.pores()
    Ts = pn.throats()
    geo = bwfpnm.Geometry.MultiPoreThroat(network=pn, name='MultiPoreThroat',
                                          pores=Ps, throats=Ts)
    bwfpnm.OpenPNM.Utilities.IO.VTK.save(pn, filename='MultiPoreThroat')
