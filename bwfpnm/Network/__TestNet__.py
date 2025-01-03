"""
===============================================================================
TestNet: Generate simple cubic network for testing purposes
===============================================================================

"""

import scipy as sp
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)

class TestNet(GenericNetwork):
    r"""
    A small cubic network for quick testing purposes

    Parameters
    ----------
        Nx, Ny, Nz: Number of nodes in x,y, and z directions
        Lc: Lattice constant -> homogeneous distance between two centers of nodes
    """

    def __init__(self, size=None, Nx=3, Ny=3, Nz=3, Lc=1, **kwargs):
        super(TestNet, self).__init__(**kwargs)
        if size is not None:
            try:
                Nx, Ny, Nz = size
            except:
                Nx, Ny, Nz = size, size, size
        self.generate(Nx, Ny, Nz, Lc)

    def generate(self, Nx, Ny, Nz, Lc):
        '''
        Create test network, of cubic geometry [5,5,5]

        Parameters
        ----------
        This network type accepts no arguments
        '''
        self._generate_setup(Nx, Ny, Nz, Lc)
        self._generate_pores()
        self._generate_throats()
        self._add_labels()
        self._macro_size()
        self._create_boundary_labels(['front', 'top', 'right', 'internal'])
        return self

    def _generate_setup(self, Nx, Ny, Nz, Lc):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._Nx = Nx
        self._Ny = Ny
        self._Nz = Nz
        self._Lc = Lc
        self._Lx = sp.float16(self._Nx*self._Lc)
        self._Ly = sp.float16(self._Ny*self._Lc)
        self._Lz = sp.float16(self._Nz*self._Lc)

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = sp.arange(0, Np)
        self['pore.all'] = sp.ones_like(ind, dtype=bool)
        pore_coords = Lc/2 + Lc*sp.array(sp.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=sp.float64).T
        self['pore.coords'] = pore_coords

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Np = Nx*Ny*Nz
        ind = sp.arange(0,Np)
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
        self['throat.all'] = sp.ones_like(sp.arange(0,sp.shape(tpore1)[0]),dtype=bool)
        self['throat.conns'] = connections

    def _add_labels(self):
        coords = self['pore.coords']
        self['pore.front'] = self.tomask(coords[:,0]<=self._Lc)
        self['pore.left'] = self.tomask(coords[:,1]<=self._Lc)
        self['pore.bottom'] = self.tomask(coords[:,2]<=self._Lc)
        self['pore.back'] = self.tomask(coords[:,0]>=(self._Lc*(self._Nx-1)))
        self['pore.right'] = self.tomask(coords[:,1]>=(self._Lc*(self._Ny-1)))
        self['pore.top'] = self.tomask(coords[:,2]>=(self._Lc*(self._Nz-1)))
        self['pore.inlet'] = self['pore.left']
        self['pore.outlet'] = self['pore.right']
        for item in ['top','bottom','left','right','front','back']:
            ps = self.pores(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self['pore.'+item][ps[:,0]]
            ps1 = self['pore.'+item][ps[:,1]]
            ts = ts[ps1*ps0]
            self['throat.'+item] = self.tomask(throats=ts)

    def domain_volume(self):
        Lc = self._Lc
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        v_mat = Nx*Lc*Ny*Lc*Nz*Lc
        return v_mat

    def _macro_size(self):
        self._macro_Lx = self._Lc
        self._macro_Ly = self._Lc
        self._macro_Lz = self._Lc

    def label_inlet_outlet(self):
        coord = self['pore.coords']
        x = coord[:, 0]
        self['pore.inlet'] = x==x.min()
        self['pore.outlet'] = x==x.max()

    def _create_boundary_labels(self,
                                labels=['front', 'top', 'right', 'internal']):
        x, y, z = self['pore.coords'].T
        xin, xout = x[self['pore.inlet']]  , x[self['pore.outlet']]
        eps = sp.amax([xin.max()-xin.min(), xout.max()-xout.min()])

        # X-direction: front-back
        if ('front' in labels) or ('back' in labels):
            try:    # for stitching purposes
                self['pore.front']
                self['pore.back']  = x <= x.min() + eps#x-rad <= x.min() #
                self['pore.front'] = x >= x.max() - eps#x+rad >= x.max() #
            except:
                self['pore.front'] = self['pore.outlet']
                self['pore.back']  = self['pore.inlet']
            self['pore.internal'] = self['pore.front'] + self['pore.back']
        else:
            self.del_properties(['pore.front', 'pore.back'])
        # Y-direction: right-left
        if ('right' in labels) or ('left' in labels):
            self['pore.left']  = y <= y.min()+eps#y-rad <= y.min()   #
            self['pore.right'] = y >= y.max()-eps#y+rad >= y.max()   #
            self['pore.internal'] += self['pore.left'] + self['pore.right']
        else:
            self.del_properties(['pore.left', 'pore.right'])
        # Z-direction: top-bottom
        if ('top' in labels) or ('bottom' in labels):
            self['pore.bottom'] = z <= z.min()+eps#z-rad <= z.min()  #
            self['pore.top']    = z >= z.max()-eps#z+rad >= z.max()  #
            self['pore.internal'] += self['pore.top'] + self['pore.bottom']
        else:
            self.del_properties(['pore.top', 'pore.bottom'])

        self['pore.internal'] = ~self['pore.internal']
        if 'internal' not in labels:
            self.del_properties(['pore.internal'])


if __name__ == '__main__':
    import bwfpnm
    pn = bwfpnm.Network.TestNet()
    print(pn.name)
