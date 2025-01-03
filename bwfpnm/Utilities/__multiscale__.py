# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 20:13:46 2016

@author: islah

===============================================================================
Network.tools.topology: Assorted topological manipulation methods
===============================================================================

"""
import scipy as _sp
#import scipy.sparse as _sprs
#import scipy.spatial as _sptl
from OpenPNM.Base import logging as _logging
from OpenPNM.Base import Controller as _controller
logger = _logging.getLogger(__name__)
_ctrl = _controller()


class multiscale(object):

    def create_macronet(self, network, min_diameter=None, shape=None):
        r'''
        Create a new network consisted of the given network's macropores
        (if min_diameter is supplied), or a network of regular lattice
        (if shape=[Nx, Ny, Nz] is supplied).
        '''


    def _voronoi_edges(self, pores):
        r'''
        Create a voronoi edges in between 2 pores as a location for subnet.
        '''

    def _subnet_positions(self, network, voronoi_edges):
        r'''
        Transform the supplied voronoi edges into cubes of subnet.
        '''

    def create_subnet(self, network):
        r'''
        Create micro/subnets from the given full network.
        '''

    def link_multinets(self, network1, network2, network3=None):
        r'''
        Create a link between networks of multiple scales.
        '''
