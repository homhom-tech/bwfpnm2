# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:15:59 2015

@author: islah

===============================================================================
Topology -- A standard 'pore & throat' geometrical model for a topological
network based on given data
===============================================================================

"""
import scipy as _sp
from bwfpnm.Geometry import models as gm
from bwfpnm.Geometry import GenericGeometry


def _porediameter(pradius, **kwargs):
    return pradius*2


def _poreproperty(pproperty, **kwargs):
    return pproperty


def _throatdiameter(tradius, **kwargs):
    return tradius*2


def _throatproperty(tproperty, **kwargs):
    return tproperty


class Topology(GenericGeometry):
    r"""
    Default geometry of Pore-throat model for topological network with given
    geometry data. The format of the data follows those of Oren&Bakke, ICL.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.
    pradius=0.0, pvolume=0.0, pshapefactor=0.0, pclayvolume=0.0,
    pconnectivity=0.0, tradius=0.0, tvolume=0.0, tshapefactor=0.0,
    tlength=0.0, tclayvolume=0.0, tporelengths=0.0, tlengthtotal=0.0

    """

    def __init__(self, pradius=1.0, pvolume=1.0, pshapefactor=1.0,
                 pclayvolume=0.0, pconnectivity=1.0,
                 tradius=1.0, tvolume=1.0, tshapefactor=1.0, tlength=1.0,
                 tclayvolume=0.0, tporelengths=1.0, tlengthtotal=1.0,
                 Gmethod='range', Amethod='G', **kwargs):
        r"""
        Initialize
        """
        super().__init__(**kwargs)
        self._generate(pradius, pvolume, pclayvolume, pshapefactor,
                       pconnectivity,
                       tradius, tvolume, tlength, tshapefactor, tporelengths,
                       tlengthtotal, tclayvolume, Gmethod, Amethod)

    def _generate(self, pradius, pvolume, pclayvolume, pshapefactor,
                  pconnectivity,
                  tradius, tvolume, tlength, tshapefactor, tporelengths,
                  tlengthtotal, tclayvolume, Gmethod, Amethod, **kwargs):
        ## Pore
        self.models.add(propname='pore.diameter',
                       model=_porediameter,
                       pradius=pradius)
        self.models.add(propname='pore.volume',
                       model=_poreproperty,
                       pproperty=pvolume)
        self.models.add(propname='pore.clayvolume',
                       model=_poreproperty,
                       pproperty=pclayvolume)
        self.models.add(propname='pore.shapefactor',
                       model=_poreproperty,
                       pproperty=pshapefactor)
        self.models.add(propname='pore.shapefactor_constant',
                        model=gm.shapefactor.constant,
                        shapefactor=self['pore.shapefactor'],
                        method=Gmethod)
        self.models.add(propname='pore.radius_eq',  # r s.t. V the same
                        model=gm.radius_eq.radius,
                        shape='arbitrary')
        self.models.add(propname='pore.area_eq',
                        model=gm.radius_eq.area)

        ## Throat
        self.models.add(propname='throat.diameter',
                       model=_throatdiameter,
                       tradius=tradius)
        self.models.add(propname='throat.length',
                       model=_throatproperty,
                       tproperty=tlength)
        self.models.add(propname='throat.volume',
                       model=_throatproperty,
                       tproperty=tvolume)
        self.models.add(propname='throat.clayvolume',
                       model=_throatproperty,
                       tproperty=tclayvolume)
        self.models.add(propname='throat.shapefactor',
                       model=_throatproperty,
                       tproperty=tshapefactor)
        self.models.add(propname='throat.porelengths',  # ordered by 'throat.conns'
                       model=_throatproperty,
                       tproperty=tporelengths)
        self.models.add(propname='throat.lengthtotal',
                       model=_throatproperty,
                       tproperty=tlengthtotal)
        if Amethod.lower() == 'g':
            self.models.add(propname='pore.area',
                           model=gm.area.shapefactor,
                           shapefactor='pore.shapefactor',
                           diameter='pore.diameter')
            self.models.add(propname='throat.area',
                           model=gm.area.shapefactor,
                           shapefactor='throat.shapefactor',
                           diameter='throat.diameter')
        elif Amethod.lower() == 'vl':
            self.models.add(propname='pore.area',
                           model=gm.area.VL,
                           volume='pore.volume',
                           length='throat.porelengths')
            self.models.add(propname='throat.area',
                           model=gm.area.VL,
                           volume='throat.volume',
                           length='throat.length')

        self.models.add(propname='throat.shapefactor_constant',
                        model=gm.shapefactor.constant,
                        shapefactor=self['throat.shapefactor'],
                        method=Gmethod)
        self.models.add(propname='throat.radius_eq',
                        model=gm.radius_eq.radius,
                        volume='throat.volume',
                        length='throat.length',
                        shapefactor='throat.shapefactor',
                        shape='arbitratry')
        self.models.add(propname='throat.area_eq',
                        model=gm.radius_eq.area,
                        volume='throat.volume',
                        length='throat.length')

    def gen_test_data(self, auto=True, prad=1.0, pvol=1.0, pclayvol=0, pG=1.0,
                       pconnectivity=0,
                       trad=1.0, tvol=1.0, tlength=1.0, tG=1.0,
                       tporelengths=1.0,
                       tlengthtotal=1.0, tclayvol=0, Gmethod='range'):
        if auto:
#            Np, Nt = self.Np, self.Nt
            prad, trad = 1.0, 1.0
            pG = 1/4/_sp.pi
            tG = 1/4/_sp.pi
            pvol, tvol = 1.0, 1.0
            pclayvol, tclayvol = 0, 0
            tlength = 1.0*_sp.ones(self.Nt)
            tporelengths = 0.1*_sp.ones((self.Nt, 2))
            tlengthtotal = tlength+tporelengths.sum(axis=1)
            pconnectivity = 0.0

        self._generate(prad, pvol, pclayvol, pG,
                       pconnectivity,
                       trad, tvol, tlength, tG, tporelengths,
                       tlengthtotal, tclayvol, Gmethod)

    def count_shape(self):
        try:
            k_Gp = self['pore.shapefactor_constant']
            k_Gt = self['throat.shapefactor_constant']
        except:
            raise Exception('shape factor data is not available')

        shape = {}
        shape['pore.Np'] = 0
        shape['throat.Nt'] = 0
        k = {'circular': 0.5, 'square': 0.5623, 'triangular': 0.6}
        for key, val in k.items():
            self['pore.'+key] = _sp.isclose(k_Gp, val)
            self['throat.'+key] = _sp.isclose(k_Gt, val)
            shape['pore.'+key] = self['pore.'+key].sum()
            shape['throat.'+key] = self['throat.'+key].sum()

            shape['pore.Np'] += shape['pore.'+key]
            shape['throat.Nt'] += shape['throat.'+key]

        # Check the total numbers of pore and throat shapes = pn.Np & Nt
        try:
            assert self._net.Np == shape['pore.Np']
        except:
            raise Exception('The number of pore shapes does not equal \
                            the number of pores')
        try:
            assert self._net.Nt == shape['throat.Nt']
        except:
            raise Exception('The number of throat shapes does not equal \
                            the number of throats')
        self._shape = shape
        return shape

    def scaling(self, scale=1, replacenet=True):
        r'''Scaling the geometry by a scaling matrix [sx, sy, sz].
        This is an in-place operation!
        NOTE: only isotropic scaling is implemented for geometric properties!
        ==> scale = a constant

        Let s = scaling constant, p & q = points [array], p' & q' = scaled p-q
        Then:
        - point: p'=sp, q'=sq
        - distance: d(p', q') = s d(p,q)
        - area: A' = s^2 A
        - volume: V' = s^3 V

        Arguments:
        ----------
        scalearray     : an array of 3 scaling components for 3 directions x, y, and z: [sx, sy, sz].
                         This array is stored in self._scalearray.
        replace         : Boolean. True -> replace the network properties with the scaled ones. False -> return (coords, [Lx, Ly, Lz])
        '''
        prop1d = ['pore.radius_eq', 'throat.radius_eq', 'throat.lengthtotal',
                  'throat.length', 'pore.diameter', 'throat.porelengths',
                  'throat.diameter']
        prop2d = ['throat.area_eq', 'throat.area', 'pore.area', 'pore.area_eq']
        prop3d = ['throat.clayvolume', 'pore.clayvolume', 'throat.volume',
                  'pore.volume']
        for prop in prop1d:
            self[prop] *= scale
        for prop in prop2d:
            self[prop] *= scale**2
        for prop in prop3d:
            self[prop] *= scale**3

        if replacenet:
            # scale the network topology, if desired
            self._net.scaling(scalearray=scale, replace=replacenet)
        self._scale = scale
        return

    def clone_geometry(self, t_stitched=None):
        r'''
        Define sizes for cloned & extended pores/throats as a result from
        Network.Topology.clone_network()
        '''
        net = self._net
        # adjust pore/throat.all labels
        del self['pore.all'], self['throat.all']
        self['pore.all'] = _sp.ones((net.Np,), dtype=bool)
        self['throat.all'] = _sp.ones((net.Nt,), dtype=bool)
        # adjust the geometry locations registered in the network
        del net['pore.'+self.name], net['throat.'+self.name]
        net['pore.'+self.name] = _sp.ones((net.Np,), dtype=bool)
        net['throat.'+self.name] = _sp.ones((net.Nt,), dtype=bool)

        # clone geometry's properties/sizes
        if t_stitched is None:
            t_stitched = net.throats('stitched')
        Nt_stitched = t_stitched.size
        properties = []
        for prop in self.props():
            prop_split = prop.split('.')
            if prop_split[0] == 'pore':
                self[prop] = _sp.hstack((self[prop], self[prop]))
            else:
                if prop_split[1] == 'porelengths':
                    continue
                # geo for the stitched pres
                if prop_split[1] in ['diameter', 'clayvolume', 'shapefactor']:
                    temp_x = _sp.array([self[prop].mean()]*Nt_stitched)
                elif prop_split[1] in ['shapefactor_constant']:
                    temp_x = gm.shapefactor.constant(self,
                                                     self['throat.shapefactor'].mean())
                    temp_x = _sp.asarray([temp_x]*Nt_stitched)
                else:
                    properties.append(prop)
                    temp_x = []
                self[prop] = _sp.hstack((self[prop], self[prop], temp_x))

        conns = net['throat.conns'][t_stitched]
        coords = net['pore.coords'][conns]
        rad = self['throat.diameter'][t_stitched]/2
        sfactor = self['throat.shapefactor'][t_stitched]

        lengthtotal = _sp.linalg.norm(_sp.diff(coords, axis=1), axis=2).flatten()
        porelengths = self['pore.diameter'][conns]/2
        length = lengthtotal - _sp.sum(porelengths, axis=1)
        area = rad**2/sfactor/4
        volume = area*length

        dict_pair = {'throat.length': length, 'throat.volume': volume,
                     'throat.lengthtotal': lengthtotal, 'throat.area': area}
        try:
            properties.remove('throat.radius_eq')
            properties.remove('throat.area_eq')
            del self['throat.area_eq'], self['throat.radius_eq']
        except:
            pass
        for prop in properties:
            self[prop] = _sp.hstack((self[prop], self[prop], dict_pair[prop]))

        prop = 'throat.porelengths'
        self[prop] = _sp.vstack((self[prop], self[prop], porelengths))

    def update_boundary_labels(self,
                                labels=['front', 'top', 'right', 'internal']):
        net = self._net
        x, y, z = net['pore.coords'].T
        rad = self['pore.diameter']/2
        eps = 3e-5  # defined by trial&error to get ~the same number as inoutlet
        try:
            eps *= self._scale
        except:
            pass

        # X-direction: front-back
        if ('front' in labels) or ('back' in labels):
            try:
                net['pore.front'] = net['pore.outlet_ori']
                net['pore.back']  = net['pore.inlet_ori']
            except:
                net['pore.front'] = net['pore.outlet']
                net['pore.back']  = net['pore.inlet']
            net['pore.internal'] = net['pore.front'] + net['pore.back']
        else:
            net.del_properties(['pore.front', 'pore.back'])
        # Y-direction: right-left
        if ('right' in labels) or ('left' in labels):
            net['pore.left']  = y-rad <= y.min()+eps   #y <= y.min()+eps
            net['pore.right'] = y+rad >= y.max()-eps   #y >= y.max()-eps
            net['pore.internal'] += net['pore.left'] + net['pore.right']
        else:
            net.del_properties(['pore.left', 'pore.right'])
        # Z-direction: top-bottom
        if ('top' in labels) or ('bottom' in labels):
            net['pore.bottom'] = z-rad <= z.min()+eps  #z <= z.min()+eps
            net['pore.top']    = z+rad >= z.max()-eps  #z >= z.max()-eps
            net['pore.internal'] += net['pore.top'] + net['pore.bottom']
        else:
            net.del_properties(['pore.top', 'pore.bottom'])

        net['pore.internal'] = ~net['pore.internal']
        if 'internal' not in labels:
            net.del_properties(['pore.internal'])

    def perturb_data(self, ratio=0.05):
        r'''
        Perturb the geometry data with uniform random error:
            [-ratio*max, ratio*max], max = max(self[property])
        '''
        for prop in self.props():
            deprop = self[prop]
            Nsize = deprop.size
            try:
                err = _sp.random.normal(0, ratio, Nsize)
            except:
                print('Property: {}'.format(prop))
                err = _sp.zeros_like(deprop)
            err = _sp.reshape(err, deprop.shape)
            deprop += err * deprop
            self[prop] = _sp.absolute(deprop)

        props = ['shapefactor_constant', 'area']
        for item in ['pore', 'throat']:
            for prop in props:
                cprop = '.'.join([item, prop])
                self.models.regenerate(cprop)


if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
