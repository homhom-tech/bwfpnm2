# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:19:52 2015

@author: islah

===============================================================================
Cylinder_Sphere -- A standard 'cylindrical & spherical' geometrical model
===============================================================================

"""
from bwfpnm.Geometry import models as bgm
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry
import scipy as sp


class Cylinder_Sphere(GenericGeometry):
    r"""
    Cylindrical and spherical subclass of GenericGeometry.
    Create a geometry objects with random diameters of pores and throats
    generated from Weibull distribution.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    scale : string
        Either 'micro' (log_range = [-9, -6]),
        'macro' (default, log_range = [-6, -3]), or
        'whole' (log_range = [-9, -3]).

    """

    def __init__(self, scale='macro', **kwargs):
        super().__init__(**kwargs)
        self._generate(scale, **kwargs)

    def _generate(self, scale, dominant_cluster=False, cluster_case=1,
                  log_range=None, **kwargs):
        if log_range is None:
            if scale == 'micro':
                log_range = [-9, -6]
            elif scale == 'whole':
                log_range = [-9, -3]
            else:
                log_range = [-6, -3]

        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        num_range=[0, 1])
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        mode='min')

        self.models.add(propname='pore.diameter',
                        model=bgm.pore_diameter.sphere2,
                        psd_name='weibull_min',
                        psd_shape=2.5,
                        psd_loc=0,
                        psd_scale=0.5,
                        log_range=log_range)
        self.models.add(propname='throat.diameter',
                        model=bgm.throat_diameter2.cylinder2,
                        tsd_name='weibull_min',
                        tsd_shape=2.5,
                        tsd_loc=0,
                        tsd_scale=0.5,
                        log_range=log_range)

        if dominant_cluster:
            self.create_dom_cluster(case=cluster_case)

        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.sphere)

        self.models.add(propname='throat.length',
                        model=gm.throat_length.straight,
                        L_negative=1e-6)   # geometric length based on coords
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.cylinder)
        self.models.add(propname='throat.area',
                        model=gm.throat_area.cylinder)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.cylinder)

        # Additional properties like in Topology network class
        ppermtr = self['pore.diameter']*sp.pi
        tpermtr = self['throat.diameter']*sp.pi
        self['pore.shapefactor'] = self['pore.area']/(ppermtr**2)
        self['throat.shapefactor'] = self['throat.area']/(tpermtr**2)
        self.models.add(propname='pore.shapefactor_constant',
                        model=bgm.shapefactor.constant,
                        shapefactor=self['pore.shapefactor'])
        self.models.add(propname='throat.shapefactor_constant',
                        model=bgm.shapefactor.constant,
                        shapefactor=self['throat.shapefactor'])
        self['throat.porelengths'] = self._calc_porelengths()

    def _calc_porelengths(self):
        conns = self['throat.conns']
        porelengths = self['pore.diameter'][conns]/2
        return porelengths

    def create_dom_cluster(self, case=1, size=None):
        r'''
        Create a dominant cluster by replacing the existing pore & throat sizes.
        Case:
        1. a straight spanning cluster
        2. a staircase spanning cluster
        '''
        td = self['throat.diameter'].max()*2
        pd = self['pore.diameter'].max()*2
        pn = self._net
        try:
            nx, ny, nz = size
        except:
            nx, ny, nz = pn._shape
        if case == 1:
            # a straight dominant spanning cluster
            pind = int(ny*sp.floor(nx/2))
            pind = sp.arange(pind, pind+ny)
            tind = int((ny-1)*sp.floor(nx/2))
            tind = sp.arange(tind, tind+ny-1)
        elif case == 2:
            # a double-straight dominant spanning cluster
            ind = int(ny*sp.floor(nx/2))
            pind = sp.arange(ind, ind + 2*ny)
            ind = int((ny-1)*sp.floor(nx/2))
            tind = sp.arange(ind, ind + 2*(ny-1))
            tind_v = pn.find_connecting_throat(pind[:ny], pind[ny:])
            tind = sp.r_[tind, sp.array(tind_v).flatten()]
            # tind must be added by the vertical edges
        elif case == 3:
            # a staircase dominant spanning cluster
            pind = sp.arange(0, nx*ny-1, ny+1)
            pind = sp.r_[pind, pind+1]
            pind.sort()
            tind = pn.find_connecting_throat(pind[:-1], pind[1:])
            tind = sp.array(tind).flatten()
        elif case == 4:
            # homogeneous size
            pind = pn.pores()
            tind = pn.throats()

        self['pore.diameter'][pind] = pd
        self['throat.diameter'][tind] = td
        self['pore.dominant_cluster'] = self.tomask(pores=pind)
        self['throat.dominant_cluster'] = self.tomask(throats=tind)


if __name__ == '__main__':
    import bwfpnm as bpnm

    ctrl = bpnm.Base.Controller()
    ctrl.loglevel = 40

    pn = bpnm.Network.RegularLattice(name='TestMicroNet', shape=[2, 2, 2])
    pn.add_boundaries(labels=['left', 'right'])

    Ps = pn.pores()
    Ts = pn.throats()
    geo = bpnm.Geometry.Cylinder_Sphere(network=pn, pores=Ps, throats=Ts,
                                        name='geo_spherical', scale='micro')
