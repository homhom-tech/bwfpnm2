# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 19:46:28 2015

@author: islah

===============================================================================
TestMultiNet: Generate multiscale cubic network for testing purposes
===============================================================================

"""

import scipy as sp
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)

class TestMultiNet(GenericNetwork):
    r"""
    A small nested multiscale cubic network for quick testing purposes

    Parameters
    ----------
        Nx, Ny, Nz: Number of nodes in x,y, and z directions
        Lc: Lattice constant -> homogeneous distance between two centers of nodes
    """

    def __init__(self, size=None, Nx=3, Ny=3, Nz=3, Lc=1,
                 micro=True, size_micro=None, fullnet=True, **kwargs):
        super().__init__(**kwargs)
        if size is not None:
            try:
                Nx, Ny, Nz = size
            except:
                Nx, Ny, Nz = size, size, size
        else:
            Nx, Ny, Nz = 2, 2, 2

        if size_micro is None:
            size_micro = [Nx, Ny, Nz]

        self.generate(Nx, Ny, Nz, Lc, size_micro)

    def generate(self, Nx, Ny, Nz, Lc, size_micro):
        '''
        Create test network, of cubic geometry [5,5,5]

        Parameters
        ----------
        This network type accepts no arguments
        '''
        self._generate_setup(Nx, Ny, Nz, Lc)
        self._micro_setup(size_micro)
        self._generate_pores()
        self._add_labels()

        return self

    def _generate_setup(self, Nx, Ny, Nz, Lc):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._Nx = Nx
        self._Ny = Ny
        self._Nz = Nz
        self._Lc = Lc
        self._Lx = sp.float16((self._Nx-1)*self._Lc)
        self._Ly = sp.float16((self._Ny-1)*self._Lc)
        self._Lz = sp.float16((self._Nz-1)*self._Lc)
        self._shape = [Nx, Ny, Nz]

    def _micro_setup(self, size_micro):
        try:
            self._micro_Nx = size_micro[0]
            self._micro_Ny = size_micro[1]
            self._micro_Nz = size_micro[2]
        except:
            self._micro_Nx = self._Nx
            self._micro_Ny = self._Ny
            self._micro_Nz = self._Nz
        # macro throat is always in micro y direction
        self._micro_Lc = self._Lc/(self._micro_Ny + 1)
        self._micro_Lx = sp.float16((self._micro_Nx-1)*self._micro_Lc)
        self._micro_Ly = sp.float16((self._micro_Ny-1)*self._micro_Lc)
        self._micro_Lz = sp.float16((self._micro_Nz-1)*self._micro_Lc)
        self._micro_shape = [self._micro_Nx, self._micro_Ny, self._micro_Nz]

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        mNx = self._micro_Nx
        mNy = self._micro_Ny
        mNz = self._micro_Nz
        mLc = self._micro_Lc
        mNp = mNx*mNy*mNz

        ind = sp.arange(0, Np)
        t_all, t_conns = self._generate_throats(Nx, Ny, Nz, Np, ind)
        Nt = sp.size(t_all)
        Np_multi = Np + Nt*mNp
        ind_multi = sp.arange(0, Np_multi)
        mind = sp.arange(0, mNp)
        self['pore.all'] = sp.ones_like(ind_multi, dtype=bool)
        self['pore.macro'] = self.tomask(sp.ones_like(ind, dtype=bool))
        self['pore.micro'] = ~self['pore.macro']

        # if one changes Lc/2 below, one have also change Lc/2 in add label function
        pore_coords = Lc/2 + Lc*sp.array(
            sp.unravel_index(ind, dims=(Nx, Ny, Nz), order='C'),
            dtype=sp.float64).T
        dims = sp.array([(mNy, mNz, mNx), (mNx, mNy, mNz), (mNz, mNx, mNy)])
        mpore_coords = []
        mt_conns = []
        for i, conns in enumerate(t_conns):
            # micro pore coordinates relative to corresponding macro pores
            coords = pore_coords[conns]     # coords of 2 connected pores
            mean_coords = sp.mean(coords, axis=0)
            direction = coords[0] != mean_coords

            mcoords = coords[0] + mLc*sp.array(
                sp.unravel_index(mind, dims=dims[direction][0], order='F'),
                dtype=sp.float64).T + mLc*direction - mLc/2*~direction
#            mcoords *= 0.25*~direction   # scaled down towards macro throat
            distance = (mcoords-coords[0])      # distance to 'throat'
            mcoords -= 0.5*distance*~direction
            mpore_coords.extend(mcoords)

            # micro throat connections
            mt_alli, mt_connsi = self._generate_throats(mNx, mNy, mNz,
                                                        mNp, mind)
            mt_connsi += Np + i*mNp

            # add boundaries: throats from boundaries to macro pores
            left = coords[0][direction] + mLc
            right = coords[0][direction] + dims[direction][0][0]*mLc
            startp = mcoords[:, direction] <= left
            stopp = mcoords[:, direction] >= right
            pnumber = sp.unique(mt_connsi)
            mpore1 = pnumber[startp[:,0]]
            mpore2 = pnumber[stopp[:,0]]
            pore1, pore2 = conns
            mt_conns_bc = [[pore1, pore] for pore in mpore1]
            mt_conns_bc.extend([[pore2, pore] for pore in mpore2])

            mt_conns.extend(mt_conns_bc)
            mt_conns.extend(mt_connsi)

        mpore_coords = sp.array(mpore_coords)
        pore_coords_all = sp.r_[pore_coords, mpore_coords]

        mt_conns = sp.array(mt_conns)
        conns = sp.r_[t_conns, mt_conns]

        self['pore.coords'] = pore_coords_all
        self['pore.micro_coords'] = pore_coords_all[self['pore.micro']]
        self['pore.macro_coords'] = pore_coords_all[self['pore.macro']]

        self['throat.all'] = sp.ones_like(sp.arange(0, sp.shape(conns)[0]),
                                          dtype=bool)
        self['throat.macro'] = self.tomask(throats=sp.ones_like(
                                           sp.arange(0, Nt), dtype=bool))
        self['throat.micro'] = ~self['throat.macro']
        self['throat.conns'] = conns


    def _generate_throats(self, Nx, Ny, Nz, Np, ind):
        r"""
        Generate the macro throats (connections, numbering and types)
        """
        #Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = sp.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = sp.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = sp.vstack((tpore1,tpore2)).T
        connections = connections[sp.lexsort((connections[:, 1], connections[:, 0]))]
        throat_all = sp.ones_like(sp.arange(0,sp.shape(tpore1)[0]),dtype=bool)
        throat_conns = connections
        return (throat_all, throat_conns)

    def _add_labels(self):
        coords = self['pore.coords'][self['pore.macro']]
        Lc = self._Lc
        self['pore.front'] = self.tomask(coords[:,0]<=Lc/2)
        self['pore.left'] = self.tomask(coords[:,1]<=Lc/2)
        self['pore.bottom'] = self.tomask(coords[:,2]<=Lc/2)
        self['pore.back'] = self.tomask(coords[:,0]>=(Lc*(self._Nx)-Lc/2))
        self['pore.right'] = self.tomask(coords[:,1]>=(Lc*(self._Ny)-Lc/2))
        self['pore.top'] = self.tomask(coords[:,2]>=(Lc*(self._Nz)-Lc/2))
        self['pore.inlet'] = self['pore.left']
        self['pore.outlet'] = self['pore.right']
        loc = ['top', 'bottom', 'left', 'right',
               'front', 'back', 'inlet', 'outlet']
        for item in loc:
            ps = self.pores(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self['pore.'+item][ps[:, 0]]
            ps1 = self['pore.'+item][ps[:, 1]]
            ts = ts[ps1*ps0]
            self['throat.'+item] = self.tomask(throats=ts)


if __name__ == '__main__':
    import bwfpnm
    pn = bwfpnm.Network.TestMultiNet(name='TestMultiNet', size=[2,2,2])
    print(pn.name)
#    pn['pore.coords'][pn['pore.macro']]
#    pn['pore.coords'][pn['pore.micro']]
#    pn['throat.conns'][pn['throat.macro']]
#    pn['throat.conns'][pn['throat.micro']]
#    bwfpnm.OpenPNM.Utilities.IO.VTK.save(pn, filename='TestMultiNet2')
