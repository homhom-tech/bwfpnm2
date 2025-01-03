# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 13:18:17 2016

@author: islah

===============================================================================
MultiNet: Generate multiscale cubic network
===============================================================================

"""

import scipy as sp
from bwfpnm.Network import TestMultiNet
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)

class MultiNet(TestMultiNet):
    r"""
    A small nested multiscale cubic network for quick testing purposes

    Parameters
    ----------
        Nx, Ny, Nz: Number of nodes in x,y, and z directions
        Lc: Lattice constant -> homogeneous distance between two centers of nodes
    """

    def __init__(self, size=None, Nx=3, Ny=3, Nz=3, Lc=1,
                 micro=True, size_micro=None, fullnet=True, **kwargs):

        super().__init__(size=size, Nx=Nx, Ny=Ny, Nz=Nz, Lc=Lc,
                         micro=micro, size_micro=size_micro, fullnet=fullnet,
                         **kwargs)


if __name__ == '__main__':
    import bwfpnm
    pn = bwfpnm.Network.MultiNet(name='MultiNet', size=[2,2,2])
    print(pn.name)
#    pn['pore.coords'][pn['pore.macro']]
#    pn['pore.coords'][pn['pore.micro']]
#    pn['throat.conns'][pn['throat.macro']]
#    pn['throat.conns'][pn['throat.micro']]
#    bwfpnm.OpenPNM.Utilities.IO.VTK.save(pn, filename='TestMultiNet2')

