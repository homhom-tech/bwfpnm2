# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 10:51:40 2017

@author: islah
===============================================================================
module __DynamicFlow__: Dynamic fluid flow
===============================================================================

Refs:
[1] C. Qin, “Water Transport in the Gas Diffusion Layer of a Polymer
Electrolyte Fuel Cell: Dynamic Pore-Network Modeling,” J. Electrochem. Soc.,
vol. 162, no. 9, pp. F1036–F1046, 2015.
[2] V. Joekar-Niasar, S. M. Hassanizadeh, and H. K. Dahle, “Non-equilibrium
effects in capillarity and interfacial area in two-phase flow: dynamic
pore-network modelling,” J. Fluid Mech., vol. 655, pp. 38–71, 2010.
"""

import scipy as sp
from bwfpnm.Algorithms import GenericMultiscaleLinearTransport
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class DynamicFlow(GenericMultiscaleLinearTransport):
    r'''
    A subclass of GenericMultiscaleLinearTransport to simulate dynamic flow.
    IMBIBITION

    Examples
    --------


    '''

    def __init__(self, network, alg_wp, phys_vapour, phys_moisture,
                 alg_dp=None,phase_wet=None, phase_dry=None, eps=1e-6,
                 **kwargs):
        r'''
        '''
        super().__init__(**kwargs)
        self._net = network
        self._phase_wet = phase_wet
        self._phase_dry = phase_dry
        self._phys_vapour = phys_vapour
        self._phys_moisture = phys_moisture
        self._alg_wp = alg_wp
        self._alg_dp = alg_dp

        logger.info('Create '+self.__class__.__name__+' Object')

        phase = self._phase_wet
        self['throat.entry_pressure'] = phase['throat.capillary_pressure']
        self['pore.entry_pressure'] = phase['pore.capillary_pressure']
        self._eps = eps

    def setup(self, conductance, quantity, super_pore_conductance):
        r"""
        This setup provides the initial data for the solver from the provided
        properties.
        It also creates the matrices A and b.
        """
        # Assigning super_pore conductance for Neumann_group BC
        if super_pore_conductance is None:
            self.super_pore_conductance = []
        else:
            self.super_pore_conductance = super_pore_conductance
        # Providing conductance values for the algorithm from the Physics name
        if sp.size(self._phase) == 1:
            self._conductance = 'throat.'+conductance.split('.')[-1]
            self._quantity = 'pore.' + self._phase.name + '_' + \
                             quantity.split('.')[-1]
            # Check health of conductance vector
            if self._phase.check_data_health(props=self._conductance).health:
                self['throat.conductance'] = self._phase[self._conductance]
            else:
                raise Exception('The provided throat conductance has problems')
        else:
            raise Exception('The linear solver accepts just one phase.')
        # Checking for the linear terms to be added to the coeff diagonal/RHS
        diag_added_data = sp.zeros(self.Np)
        RHS_added_data = sp.zeros(self.Np)
        for label in self.labels():
            if 'pore.source_' in label:
                source_name = 'pore.' + \
                              (label.split('.')[-1]).replace('source_', '')
                matching_physics = [phys for phys in self._phase._physics
                                    if source_name in phys.models.keys()]
                for phys in matching_physics:
                    x = phys.models[source_name]['x']
                    if x != '' and type(x) == str:
                        if x.split('.')[-1] != quantity.split('.')[-1]:
                            raise Exception('The quantity(pore.' +
                                            x.split('.')[-1] +
                                            '), provided by source term(' +
                                            source_name + '), is different ' +
                                            'from the main quantity(pore.' +
                                            quantity.split('.')[-1] + ') in ' +
                                            self.name + ' algorithm.')
                source_name = label.replace('pore.source_', '')
                if 'pore.source_linear_s1_' + source_name in self.props():
                    prop1 = 'pore.source_linear_s1_' + source_name
                    pores = ~sp.isnan(self[prop1])
                    diag_added_data[pores] = diag_added_data[pores] + \
                        self[prop1][pores]
                    prop2 = 'pore.source_linear_s2_' + source_name
                    pores = ~sp.isnan(self[prop2])
                    RHS_added_data[pores] = RHS_added_data[pores] + \
                        self[prop2][pores]
        # Creating A and b based on the conductance values and new linear terms
        logger.info('Creating Coefficient matrix for the algorithm')
        d = diag_added_data
        self.A = self._build_coefficient_matrix(modified_diag_pores=self.Ps,
                                                diag_added_data=d)
        logger.info('Creating RHS matrix for the algorithm')
        self.b = self._build_RHS_matrix(modified_RHS_pores=self.Ps,
                                        RHS_added_data=-RHS_added_data)

    def _build_coefficient_matrix(self, modified_diag_pores=None,
                                  diag_added_data=None, mode='overwrite'):
        r"""
        This builds the sparse coefficient matrix for the linear solver.
        """
        if mode == 'overwrite':
            # Filling coefficient matrix
            tpore1 = self._net['throat.conns'][:, 0]
            tpore2 = self._net['throat.conns'][:, 1]
            # Identify Dirichlet pores
            try:
                temp = self.pores(self._phase.name + '_Dirichlet',
                                  mode='difference')
            except KeyError:
                temp = self.Ps
                logger.warning('No direct Dirichlet boundary condition has ' +
                               'been applied to the phase ' +
                               self._phase.name + ' in the algorithm ' +
                               self.name)
            loc1 = sp.in1d(tpore1, temp)
            loc2 = sp.in1d(tpore2, temp)
            modified_tpore1 = tpore1[loc1]
            modified_tpore2 = tpore2[loc1]
            row = modified_tpore1
            col = modified_tpore2
            # Expand the conductance to a vector if necessary
            g = self['throat.conductance']
            if sp.size(g) == 1:
                g = g * sp.ones(self.Nt)
            data_main = g
            data = data_main[loc1]
            modified_tpore2 = tpore2[loc2]
            modified_tpore1 = tpore1[loc2]
            row = sp.append(row, modified_tpore2)
            col = sp.append(col, modified_tpore1)
            data = sp.append(data, data_main[loc2])
            A_dim = self.Np
            # Check for Neuman_group BCs and add superpores if necessary
            try:
                self.pores(self._phase.name + '_Neumann_group')
                self._extra_Neumann_size = len(getattr(self, '_pore_' +
                                                       self._phase.name +
                                                       '_Neumann_group_' +
                                                       'location'))
                self._group_Neumann_vals = sp.zeros(self._extra_Neumann_size)
                l_g_super = len(self.super_pore_conductance)
                if l_g_super not in [0, 1, self._extra_Neumann_size]:
                    raise Exception('length of the list of super_pore_'
                                    'conductance and the number of different'
                                    ' Neumann_group BCs do not match.')
                if l_g_super == 1:
                    t = [sp.array(self.super_pore_conductance)]
                    self.super_pore_conductance = t * self._extra_Neumann_size
                for N in sp.arange(0, self._extra_Neumann_size):
                    neu_tpore2 = getattr(self, '_pore_' + self._phase.name +
                                         '_Neumann_group_location')[N]
                    Nval = self['pore.' + self._phase.name +
                                '_bcval_Neumann_group']
                    self._group_Neumann_vals[N] = sp.unique(Nval[neu_tpore2])
                    nt = self._net.find_neighbor_throats(pores=neu_tpore2)
                    try:
                        g_super = self.super_pore_conductance[N]
                    except IndexError:
                        g_super = 1e-3 * min(data_main[nt])
                        self.super_pore_conductance.append(g_super)
                    if sp.size(g_super) == 1:
                        g_super = len(neu_tpore2)*[g_super]
                    row = sp.append(row, neu_tpore2)
                    col = sp.append(col, len(neu_tpore2) * [A_dim + N])
                    data = sp.append(data, g_super)
                    row = sp.append(row, len(neu_tpore2) * [A_dim + N])
                    col = sp.append(col, neu_tpore2)
                    data = sp.append(data, g_super)
                A_dim = A_dim + self._extra_Neumann_size
            except KeyError:
                pass
            # Adding positions for diagonal
            diag = sp.arange(0, A_dim)
            try:
                pores = self.pores(self._phase.name + '_Dirichlet')
                row = sp.append(row, diag[pores])
                col = sp.append(col, diag[pores])
                data = sp.append(data, sp.ones_like(diag[pores]))
                temp_data = sp.copy(data)
                temp_data[sp.in1d(row, diag[pores])] = 0
                non_Dir_diag = diag[~sp.in1d(diag, diag[pores])]
            except KeyError:
                temp_data = sp.copy(data)
                non_Dir_diag = diag
            S_temp = sp.zeros(A_dim)
            for i in sp.arange(0, len(row)):
                S_temp[row[i]] = S_temp[row[i]] - temp_data[i]
            # Store values for modifying the diagonal in mode='modify_diagonal'
            self._non_source_row = row
            self._non_source_col = col
            self._non_source_data = data
            self._non_Dir_diag = non_Dir_diag
            self._diagonal_vals = S_temp
            self._coeff_dimension = A_dim

        if mode in ['overwrite', 'modify_diagonal']:
            diagonal_vals = sp.copy(self._diagonal_vals)
            # Adding necessary terms to the diagonal such as source terms
            if modified_diag_pores is not None and diag_added_data is not None:
                if sp.size(modified_diag_pores) == sp.size(diag_added_data):
                    sec1 = self._diagonal_vals[modified_diag_pores]
                    sec2 = diag_added_data
                    diagonal_vals[modified_diag_pores] = sec1 + sec2
                else:
                    raise Exception('Provided data and pores for modifying '
                                    'coefficient matrix should have the same' +
                                    ' size!')
                if mode == 'overwrite':
                    self._diagonal_vals = diagonal_vals
            data = sp.append(self._non_source_data,
                             diagonal_vals[self._non_Dir_diag])
            row = sp.append(self._non_source_row, self._non_Dir_diag)
            col = sp.append(self._non_source_col, self._non_Dir_diag)
            # Convert the lists to the sparse matrix
            a = sprs.coo.coo_matrix((data, (row, col)),
                                    (self._coeff_dimension,
                                     self._coeff_dimension))
            A = a.tocsr()
            A.eliminate_zeros()
            return(A)

    def _build_RHS_matrix(self, modified_RHS_pores=None, RHS_added_data=None,
                          mode='overwrite'):
        r"""
        This builds the right-hand-side matrix for the linear solver.
        """

        if mode == 'overwrite':
            A_dim = self._coeff_dimension
            b = sp.zeros([A_dim, 1])
            try:
                Dir_pores = self.pores(self._phase.name + '_Dirichlet')
                Dir_pores_vals = self['pore.' + self._phase.name +
                                      '_bcval_Dirichlet'][Dir_pores]
                b[Dir_pores] = sp.reshape(Dir_pores_vals, [len(Dir_pores), 1])
            except KeyError:
                pass
            try:
                ind_Neu_pores = self.pores(self._phase.name + '_Neumann')
                ind_Neu_pores_vals = self['pore.' + self._phase.name +
                                          '_bcval_' + 'Neumann'][ind_Neu_pores]
                b[ind_Neu_pores] = sp.reshape(ind_Neu_pores_vals,
                                              [len(ind_Neu_pores), 1])
            except KeyError:
                pass
            try:
                self.pores(self._phase.name + '_Neumann_group')
                pnum = self._net.Np
                NG_loc = sp.r_[pnum: (pnum + len(self._group_Neumann_vals))]
                NG_l = len(self._group_Neumann_vals)
                NG_arr = self._group_Neumann_vals[sp.r_[0:NG_l]]
                b[NG_loc] = sp.reshape(NG_arr, [NG_l, 1])
            except KeyError:
                pass
        if mode == 'modify_RHS':
            b = sp.copy(self.b)
        if mode in ['overwrite', 'modify_RHS']:
            # Adding necessary terms to the RHS for non-Dirichlet pores
            if modified_RHS_pores is not None and RHS_added_data is not None:
                if sp.size(modified_RHS_pores) == sp.size(RHS_added_data):
                    p = sp.in1d(modified_RHS_pores, self._non_Dir_diag)
                    data = RHS_added_data[p]
                    b[modified_RHS_pores[p]] = b[modified_RHS_pores[p]] + \
                                               data.reshape([len(data), 1])
                else:
                    raise Exception('Provided data and pores for modifying'
                                    ' RHS matrix should have the same size!')
        return(b)