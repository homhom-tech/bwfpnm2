# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 12:40:47 2015

@author: islah

===============================================================================
module __WaterFlow__: Viscous fluid flow -> liquid permeability
===============================================================================

"""
import scipy as sp
from bwfpnm.Algorithms import GenericMultiscaleLinearTransport
from OpenPNM.Base import logging
import re
logger = logging.getLogger(__name__)


class WaterFlow(GenericMultiscaleLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the network.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,
    ...                                     pores=pn.pores(),
    ...                                     throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,
    ...                                     pores=pn.pores(),throats=pn.throats())
    >>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=phase1)
    >>> BC1_pores = pn.pores('top')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    >>> BC2_pores = pn.pores('bottom')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    >>> alg.run()
    >>> alg.return_results()
    >>> Peff = round(alg.calc_eff_permeability(), 10)
    >>> print(Peff)
    1.8663e-05
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')

    def setup(self, conductance='hydraulic_conductance', quantity='pressure',
              super_pore_conductance=None, **params):
        r"""
        This setup provides the initial requirements for the solver setup.
        """
        logger.info('Setup ' + self.__class__.__name__)
        super().setup(conductance=conductance, quantity=quantity,
                      super_pore_conductance=super_pore_conductance)

    def calc_eff_permeability(self, conductance=None):
        r"""
        This calculates the effective permeability in this linear
        transport algorithm.
        """
        self._eff_permeability = self._calc_eff_prop(conductance=conductance)
        return self._eff_permeability

    def _calc_eff_prop(self, conductance, check_health=False):
        r'''
        This returns the main parameters for calculating the effective property in a linear transport equation.
        It also checks for the proper boundary conditions, inlets and outlets.

        Parameters
        ----------
        check_health : boolean(optional)
            It analyzes the inlet and outlet pores to check their spatial positions
        '''
        try:
            self[self._quantity]
        except:
            raise Exception('The algorithm has not been run yet. Cannot calculate effective property.')
        #Determine boundary conditions by analyzing algorithm object
        try:
            Ps = self.pores('pore.'+self._phase.name+'_Dirichlet')
            BCs = sp.unique(self['pore.'+self._phase.name+'_bcval_Dirichlet'][Ps])
            if sp.shape(BCs)[0] != 2:
                raise Exception('The supplied algorithm did not have appropriate BCs')
            inlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amax(BCs))[0]
            outlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amin(BCs))[0]
        except:
            # for OpenPNM 1.4.6
            Ps = self.pores('pore.Dirichlet')
            BCs = sp.unique(self['pore.bcval_Dirichlet'][Ps])
            if sp.shape(BCs)[0] != 2:
                raise Exception('The supplied algorithm did not have appropriate BCs')
            inlets = sp.where(self['pore.bcval_Dirichlet']==sp.amax(BCs))[0]
            outlets = sp.where(self['pore.bcval_Dirichlet']==sp.amin(BCs))[0]

        #Analyze input and output pores
        if check_health:
            #Check for coplanarity
            if self._net.iscoplanar(inlets) == False:
                raise Exception('The inlet pores do not define a plane.\
                    Effective property will be approximation')
            if self._net.iscoplanar(outlets) == False:
                raise Exception('The outlet pores do not define a plane.\
                    Effective property will be approximation')
            #Ensure pores are on a face of domain (only 1 non-self neighbor each)
            PnI = self._net.find_neighbor_pores(pores=inlets,
                                                mode='not_intersection',
                                                excl_self=True)
            if sp.shape(PnI) != sp.shape(inlets):
                logger.warning('The inlet pores have too many neighbors.\
                    Internal pores appear to be selected.')
                pass
            PnO = self._net.find_neighbor_pores(pores=outlets,
                                                mode='not_intersection',
                                                excl_self=True)
            if sp.shape(PnO) != sp.shape(outlets):
                logger.warning('The outlet pores have too many neighbors.\
                    Internal pores appear to be selected.')
                pass

        #Fetch area and length of domain
        try:
            A = self._net.domain_area(face='xy')
            L = self._net.domain_length(direction='x')
        except:
            A, L = 1, 1
            logger.warning('The macroscopic length and area are set to unity')

        self._inlets = inlets
        self._outlets = outlets
        if inlets.size < outlets.size:
            plane = 'inlet'
            xplane = 'outlet'
        else:
            plane = 'outlet'
            xplane = 'inlet'
#        plane = 'inlet'
        flow = self.rate(conductance=conductance, plane=plane)
        flow2 = self.rate(conductance=conductance, plane=xplane)

#        flow_out = self.rate(pores=outlets, conductance=conductance)
#        flow = -flow_out
#        flow = flow_in
        self._BCplane = plane
        self._flow2 = flow2
        self._flow = flow

#        flow = sp.mean([flow_in, flow_out])
#        print('Boundary plane: \t' + plane)
#        print('FLow \t {}\n'.format(flow1[0]))
#        print('FLow out\t {}\n\n'.format(flow_out))
        K1 = sp.sum(flow)*L/A/(BCs[0]-BCs[1])
        K2 = sp.sum(flow2)*L/A/(BCs[0]-BCs[1])

#        self._net._permeability_xplane = K2
        return K1, K2

    def calc_abs_permeability(self):
        r'''
        Calculate the absolute intrinsic permeability [m^2] calculated from
        wet (& dry) effective permeability (K_s=1), (& K_s=0).
        The result is stored in self._abs_m2_permeability and/or self._abs_m2_permeability_dry
        '''
        phases = self._net._phases
        for phase in phases:
            words = phase._name.lower().split('_')
            words += re.findall('[A-Z][^A-Z]*', phase._name)
            words = [word.lower() for word in words]
            if 'water' in words:
                water = phase
            elif 'vapour' in words:
                vapour = phase
            elif 'vap' in words:
                vapour = phase

        # K_abs from wet state --> water phase
        try:
            mu = sp.mean(water['pore.viscosity'])
            rho = sp.mean(water['pore.density'])
            self._abs_m2_permeability = sp.amax(self._eff_permeability)*mu/rho
        except:
            logger.error('the effective permeability has not been calculated')

        # K_abs from dry state --> vapour phase
        try:
            mu = sp.mean(vapour['pore.viscosity'])
            rho = sp.mean(vapour['pore.density'])
            self._abs_m2_permeability_dry = sp.amin(self._eff_permeability)*mu/rho
        except:
            logger.warning('Failed to calculate abs permeability from\
                            dry state')
        return self._abs_m2_permeability

    def calc_mD_permeability(self):
        r'''
        Convert absolute permeability from m^2 to mD unit, and store it in
        self._abs_mD_permeability

        Requirement
        -----------
        absolute permeability or at least effective permeability has been calculated
        '''
        cons = 1e3/9.869233e-13
        try:
            # from m^2 to mD
            self._abs_mD_permeability = self._abs_m2_permeability*cons
        except:
            self._abs_mD_permeability = self.calc_abs_permeability*cons

        return self._abs_mD_permeability


    def store_result(self, Pc=None, sat=None, sat_surf=None, sat_moist=None,
                     w_sat=None, k=None, g=None, span_water=[], span_air=[]):
        r"""
        Storing Pc, sat, w and k in the object
        """
        self._Pc = sp.array(Pc)
        self._log_Pc = sp.log10(-self._Pc)
        self._saturation = sp.array(sat)
        self._moisturecontent = self._saturation*w_sat
        self._permeability = sp.array(k)
        self._conductance = sp.array(g)
        self._span_water = sp.array(span_water)
        self._span_air = sp.array(span_air)
        try:
            self._log_permeability = sp.log10(self._permeability)
        except:
            self._log_permeability = sp.zeros_like(self._Pc)

        try:
            self._saturation_surf = sp.array(sat_surf)
        except:
            self._saturation_surf = sp.zeros_like(sat)

        try:
            self._saturation_vapour = sp.array(sat_moist)
        except:
            self._saturation_vapour = sp.zeros_like(sat)
#        # save to the network instance/object
#        self._net._Pc = self._Pc
#        self._net._saturation = self._net._saturation
#        self._net._moisturecontent = self._net._moisturecontent
#        self._net._permeability = self._net._permeability

    def _modify_system(self, row_scaling=False, **kwargs):
        r"""
        self.setup() create a linear system Ax=b where:
        A[Np, Np]: a full-rank conduction matrix
        x[Np, 1]: a pressure array of all nodes, including the BCs
        b[Np, 1]: a zero rhs array except at boundary nodes.

        This function modify that system to A'x'=b', where:
        A'[N, N]: a (low-rank) conduction matrix excluding the BCs
        x'[N, 1]: a pressure array of internal nodes
        b'[N, 1]: a zero rhs array except at nodes next to Bc nodes
        N = Np - N_bc: number of internal nodes which pressures are unknown
        """
        # Identify Dirichlet pores
        try:
            inpores = self.pores(self._phase.name + '_Dirichlet',
                              mode='difference')
            bcpores = self.pores(self._phase.name + '_Dirichlet')
        except KeyError:
            inpores = self.Ps
            bcpores = []
            logger.warning('No direct Dirichlet boundary condition has ' +
                           'been applied to the phase ' +
                           self._phase.name + ' in the algorithm ' +
                           self.name)
        self._inpores = inpores
        self._bcpores = bcpores

        A = self.A
        b = self.b

        # Modify rhs b
        bb = A*b
        self.b = -bb[inpores]
        # Cut BC rows & cols of A
        self._set_row(A, bcpores, 0)
        self._set_col(A, bcpores, 0)
        A.eliminate_zeros()
        self._reshape_matrix(A, bcpores, inpores, **kwargs)
        if row_scaling:
            self._row_scaling()

    def _set_row(self, A, ind, value=0):
        r"""
        Setting rows ind of a csr matrix to the given val.
        """
#        A = self.A
        if not isinstance(A, sp.sparse.csr_matrix):
            raise ValueError('Matrix given must be of CSR format.')
        for i in ind:
            A.data[A.indptr[i]:A.indptr[i+1]] = value

    def _set_col(self, A, ind, value=0):
        r"""
        Setting columns ind of a csr matrix to the given val.
        """
#        A = self.A
        if not isinstance(A, sp.sparse.csr_matrix):
            raise ValueError('Matrix given must be of CSR format.')
        loc = sp.in1d(A.indices, ind)
        A.data[loc] = value

    def _reshape_matrix(self, A, bc, inpore,
                        indices=None, indptr=None, shape=None, **kwargs):
        r"""
        Reshape matrix after eliminating zeros in ind rows and cols
        This operation is only need to be performed once,
        since the network is the same.
        """
        if indices is not None:
            A.indices = indices
            A.indptr = indptr
            A._shape = shape
        else:
            # modify A.indices
            cols = A.indices
            temps = bc - sp.arange(bc.size)
            for temp in temps:      # this loop is slow!
                loc = cols > temp
                cols[loc] -= 1
            # modify A.indptr & A._shape
            A.indptr = sp.delete(A.indptr, bc)
            N = inpore.size
            A._shape = (N,N)

    def _row_scaling(self, **kwargs):
        r'''
        Modify matrix A by its diagonal, and rhs as well.
        '''
        v = 1./self.A.diagonal()
        nv = len(v)
        self.A = sp.sparse.spdiags(v, 0, nv, nv) * self.A

        b = self.b.T * v
        self.b = sp.reshape(b, self.b.shape)

    def solve(self, A=None, b=None, iterative_solver=None, **kwargs):
        r"""
        Executes the right algorithm for the solution: regular solution of a
        linear system or iterative solution over the nonlinear source terms.

        Parameters
        ----------
        A : sparse matrix
            2D Coefficient matrix
        b : dense matrix
            1D RHS vector
        iterative_sovler : string
            Name of solver to use.  If not solve is specified, sp.solve is used
            which is a direct solver (SuperLU on default Scipy installation)
        kwargs : list of keyword arguments
            These arguments and values are sent to the sparse solver, so read
            the specific documentation for the solver chosen
        """
        self._iterative_solver = iterative_solver

        # Executes the right algorithm
        if any('pore.source_nonlinear' in s for s in self.props()):
            X = self._do_one_outer_iteration(**kwargs)
        else:
            X = self._do_one_inner_iteration(A, b, **kwargs)
        self.X = X
        self._Neumann_super_X = self.X[self.Np:self._coeff_dimension]
        # Removing the additional super pore variables from the results
        try:
            Dir_pores_vals = self['pore.' + self._phase.name +
                                          '_bcval_Dirichlet'][self._bcpores]
            self[self._quantity] = 0.0
            self[self._quantity][self._inpores] = self.X
            self[self._quantity][self._bcpores] = Dir_pores_vals
        except:
            self[self._quantity] = self.X[self.Ps]
        logger.info('Writing the results to ' + '[\'' + self._quantity +
                    '\'] in the ' + self.name + ' algorithm.')

    def rate(self, pores=None, network=None, conductance=None, X_value=None,
             mode='group', plane=None):
        r'''
        Send a list of pores and receive the net rate
        of material moving into them.

        Parameters
        ----------
        pores : array_like
            The pores where the net rate will be calculated
        network : OpenPNM Network Object
            The network object to which this algorithm will apply.
            If no network is sent, the rate will apply to the network which is attached to the algorithm.
        conductance : array_like
            The conductance which this algorithm will use to calculate the rate.
            If no conductance is sent, the rate will use the 'throat.conductance' which is attached to the algorithm.
        X_value : array_like
            The values of the quantity (temperature, mole_fraction, voltage, ...), which this algorithm will use to calculate the rate.
            If no X_value is sent, the rate will look at the '_quantity', which is attached to the algorithm.
        mode : string, optional
            Controls how to return the rate.  Options are:
            - 'group'(default): It returns the cumulative rate moving into them
            - 'single': It calculates the rate for each pore individually.

        '''
        if network is None: network = self._net
        if conductance is None: conductance = self['throat.conductance']
        if X_value is None: X_value = self[self._quantity]
        if plane == 'inlet':
            pores = self._inlets
        elif plane == 'outlet':
            pores = self._outlets

        pores = sp.array(pores,ndmin=1)
        R = []
        if mode=='group':
            t = network.find_neighbor_throats(pores, flatten=True,
                                              mode='not_intersection')
            throat_group_num = 1
        elif mode=='single':
            t = network.find_neighbor_throats(pores, flatten=False,
                                              mode='not_intersection')
            throat_group_num = sp.size(t)

        for i in sp.r_[0:throat_group_num]:
            if mode=='group':
                throats = t
                P = pores
            elif mode=='single':
                throats = t[i]
                P = pores[i]
            p1 = network.find_connected_pores(throats)[:,0]
            p2 = network.find_connected_pores(throats)[:,1]
            pores1 = sp.copy(p1)
            pores2 = sp.copy(p2)
            #Changes to pores1 and pores2 to make them as the inner and outer pores
            pores1[-sp.in1d(p1,P)] = p2[-sp.in1d(p1,P)] #pores1 = BC pores
            pores2[-sp.in1d(p1,P)] = p1[-sp.in1d(p1,P)] #pores2 = inner pores
            X1 = X_value[pores1] # the BC Pc, a unique Pc
            X2 = X_value[pores2] # should be all smaller or bigger than X1
            g = conductance[throats]
            dPc = X2-X1
            rate = sp.multiply(g,dPc)
##            tag = dPc!=0
##            tag2 = sp.absolute(dPc) > 1e-5
            if plane == 'inlet':
                # for inlet: X1 is the max
                t_ok = X1>X2
                t_eq = X1==X2
                t_err = X1<X2
                id_ok = throats[t_ok]
                id_eq = throats[t_eq]
                id_err = throats[t_err]
            elif plane == 'outlet':
                # for inlet: X1 is the max
                t_ok = X1<X2
                t_eq = X1==X2
                t_err = X1>X2
                id_ok = throats[t_ok]
                id_eq = throats[t_eq]
                id_err = throats[t_err]
                rate *= -1
#==============================================================================
#             errAx = self.A.dot(self.X.reshape(self.b.shape)) - self.b
#             p_err = sp.where(errAx)[0]
#             X_val_err = errAx[pores2[t_err]]
#             residual2 = sp.linalg.norm(errAx, ord=2)
#             residual1 = sp.linalg.norm(errAx, ord=1)
#             relres = sp.linalg.norm(residual2/sp.absolute(self.X))
#             print('\n Bc Pc: \t{}'.format(X1[0]))
#             print('L2 residual: \t{}'.format(residual2))
#             print('L1 residual: \t{}'.format(residual1))
#             print('L2 relatif residual: \t{}'.format(relres))
#             print('Number of nonzero residual: \t{}'.format(p_err.size))
#             print('Number of error press diff: \t{}'.format(t_err.sum()))
#==============================================================================
#            print('Ok at throats IDs: \n{}\n'.format(throats[t_ok]))
#            print('Equal at throats IDs: \n{}\n'.format(throats[t_eq]))
#            print('Error at throats IDs: \n{}\n'.format(throats[t_err]))

            R.append(sp.sum(rate[t_ok]))
        return(sp.array(R, ndmin=1))

    def return_rate(self, pores=None, throats=None, case='', **kwargs):
        r"""
        Send rate of simulation to phase.
        """
        if pores is None:
            pores = self.Ps
        if throats is None:
            throats = self.Ts

        phase_quantity = self._quantity.replace(self._phase.name + '_', '')
        if phase_quantity not in self._phase.props():
            self._phase[phase_quantity] = sp.nan
        self._phase[phase_quantity][pores] = self[self._quantity][pores]
        conn_arr = self._net.find_connected_pores(self.Ts)
        dx = sp.squeeze(sp.diff(self[self._quantity][conn_arr], n=1, axis=1))
        g = self['throat.conductance']
        rate = sp.absolute(g * dx)
        if 'throat.rate' not in self._phase.props():
            self._phase['throat.rate'+case] = sp.nan
        self._phase['throat.rate'+case][throats] = rate[throats]
        self._phase['throat.delta_pressure'+case] = dx
        logger.debug('Results of ' + self.name +
                     ' algorithm have been added to ' + self._phase.name)
