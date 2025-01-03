# -*- coding: utf-8 -*-
"""
===============================================================================
module __StokesFlow__: Viscous fluid flow
===============================================================================

"""
import scipy as sp
from bwfpnm.Algorithms import GenericMultiscaleLinearTransport
from OpenPNM.Base import logging
from bwfpnm.Phases.models import vapour_pressure as pv
import re
#import numexpr as ne
#from numba import jit, vectorize, guvectorize, int32, boolean
logger = logging.getLogger(__name__)


#==============================================================================
# #%% numba - vectorize
# #==========================================
# @vectorize([boolean(int32, int32)], nopython=True, target='parallel')
# def ff(cols, temp):
#     '''a scalar operation, will be vectorized by numba'''
#     return cols > temp
#
# @vectorize([int32(int32, boolean)], nopython=True, target='parallel')
# def gf(cols, loc):
#     return cols - loc*1
#
# @jit
# def _reshape_matrix(A, bc, inpore):
#     r"""
#     Reshape matrix after eliminating zeros in ind rows and cols
#     This operation is only need to be performed once,
#     since the network is the same.
#     """
#     # modify A.indices
#     cols = A.indices
#     temps = bc - sp.arange(bc.size)
#     for temp in temps:      # this loop is slow! improved!
#         loc = ff(cols, temp)
#         cols = gf(cols, loc)
#     # modify A.indptr & A._shape
#     A.indices = cols
#     A.indptr = sp.delete(A.indptr, bc)
#     N = inpore.size
#     A._shape = (N,N)
#==============================================================================


#==============================================================================
# #%% numba - guvectorize
# #==========================================
# @jit
# def _reshape_matrix(A, bc, inpore):
#     r"""
#     Reshape matrix after eliminating zeros in ind rows and cols
#     This operation is only need to be performed once,
#     since the network is the same.
#     """
#     A.indptr = sp.delete(A.indptr, bc)
#     N = inpore.size
#     A._shape = (N,N)
#     # modify A.indices
#     temps = bc - sp.arange(bc.size)
#     cols = _modify_indices(A.indices, temps.astype('int32'))
# #    A.indices = cols
# @guvectorize([(int32[:], int32[:], int32[:])], '(n),(m)->(n)',
#                nopython=True, target='parallel')
# def _modify_indices(cols, temps, output):
#     for temp in temps:
#         i = 0
#         for col in cols:
#             if col > temp:
#                 cols[i] = col - 1
#             i += 1
#     output = cols
# #    return output
#
#==============================================================================

#%% numba - jit
#==========================================
#@jit
#def _reshape_matrix(A, bc, inpore):
#    r"""
#    Reshape matrix after eliminating zeros in ind rows and cols
#    This operation is only need to be performed once,
#    since the network is the same.
#    """
#    A.indptr = sp.delete(A.indptr, bc)
#    N = inpore.size
#    A._shape = (N,N)
#    # modify A.indices
#    temps = bc - sp.arange(bc.size)
#    _modify_indices(A, temps)
#
#@jit
#def _modify_indices(A, temps):
#    cols = A.indices
#    for temp in temps:
#        for i, col in enumerate(cols):
#            if col > temp:
#                cols[i] -= 1
#==========================================




class MoistureFlow(GenericMultiscaleLinearTransport):
    r'''
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the network.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
    >>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=phase1)
    >>> BC1_pores = pn.pores('top')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=BC1_pores)
    >>> BC2_pores = pn.pores('bottom')
    >>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.4, pores=BC2_pores)
    >>> alg.run()
    >>> alg.return_results()
    >>> Peff = round(alg.calc_eff_permeability(), 10)
    >>> print(Peff) #unless something changed with our test objects, this should print "1.8663e-05"
    1.8663e-05

    '''

    def __init__(self,**kwargs):
        r'''
        '''
        super(MoistureFlow, self).__init__(**kwargs)
        logger.info('Create '+self.__class__.__name__+' Object')

    def setup(self, conductance='conduit_conductance',
              quantity='pressure', super_pore_conductance=None, **params):
        r'''
        This setup provides the initial requirements for the solver setup.
        '''
        logger.info("Setup "+self.__class__.__name__)
        super(MoistureFlow, self).setup(
            conductance=conductance, quantity=quantity,
            super_pore_conductance=super_pore_conductance)

    def calc_eff_permeability(self, conductance=None):
        r'''
        This calculates the effective permeability [s] = [kg/msPa]
        in this linear transport algorithm.
        '''
        self._eff_permeability = self._calc_eff_prop(conductance=conductance)
        return self._eff_permeability

    def calc_eff_conduct_conceptual(self, conductance=None):
        r'''
        This calculates the effective flow rate/BC [kg/s/Pa] = [sm] conductance
        in this linear transport algorithm.
        '''
        G = self._calc_eff_prop_conceptual(conductance=conductance)
        self._eff_conductance = G
        return self._eff_conductance

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
#        try:
#            A = self._net.domain_area(face='xy')
#            L = self._net.domain_length(direction='x')
#        except:
        A, L = self._calc_domain_dimension()

        self._inlets = inlets
        self._outlets = outlets
#        if inlets.size < outlets.size:
        plane = 'inlet'
        xplane = 'outlet'
#        else:
#            plane = 'outlet'
#            xplane = 'inlet'
#        plane = 'inlet'

        flow = self.rate(conductance=conductance, plane=plane)   # inlet
        flow2 = self.rate(conductance=conductance, plane=xplane) # outlet
        print(flow)
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

        for item in self._net._phases:
            if item.name == 'vapour':
                phase = item
                break
        pv1 = pv.pore(phase, BCs[0])[0]
        pv2 = pv.pore(phase, BCs[1])[0]
        deltav1 = sp.sum(flow)*L/A/(pv1 - pv2)
        deltav2 = sp.sum(flow2)*L/A/(pv1 - pv2)
        # print(A, L)
        # print(deltav1, deltav2)
#        self._net._permeability_xplane = K2
        return K1, K2, deltav1, deltav2

    def _calc_domain_dimension(self):
        pn = self._net
        coords = pn['pore.coords']

        pinlet = pn.pores('inlet')
        poutlet = pn.pores('outlet')
        pos_inlet = sp.mean(coords[pinlet], axis=0)
        pos_outlet = sp.mean(coords[poutlet], axis=0)
        dpos = sp.absolute(pos_outlet - pos_inlet)

        direction = sp.where(dpos==dpos.max())[0]

        try:
            Lx, Ly, Lz = pn._Lx, pn._Ly, pn._Lz
        except:
            Lx, Ly, Lz = pn._macro_Lx, pn._macro_Ly, pn._macro_Lz

        try:
            if direction==1:
                A = Lx *Lz
                L = Ly
            elif direction==0:
                A = Ly * Lz
                L = Lx
            else:
                A = Lx * Ly
                L = Lz
        except:
            A, L = 1, 1
            logger.warning('The macroscopic length and area are set to unity')
        return A, L

    def _calc_eff_prop_conceptual(self, conductance, check_health=False):
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
        Ps = self.pores('pore.'+self._phase.name+'_Dirichlet')
        BCs = sp.unique(self['pore.'+self._phase.name+'_bcval_Dirichlet'][Ps])
        if sp.shape(BCs)[0] != 2:
            raise Exception('The supplied algorithm did not have appropriate BCs')
        inlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amax(BCs))[0]
        outlets = sp.where(self['pore.'+self._phase.name+'_bcval_Dirichlet']==sp.amin(BCs))[0]

        #Analyze input and output pores
        if check_health:
            #Check for coplanarity
            if self._net.iscoplanar(inlets) == False:
                raise Exception('The inlet pores do not define a plane. \
                    Effective property will be approximation')
            if self._net.iscoplanar(outlets) == False:
                raise Exception('The outlet pores do not define a plane.\
                    Effective property will be approximation')
            # Ensure pores are on a face of domain (only 1 non-self neighbor each)
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

        flow = self.rate(pores=inlets, conductance=conductance)

        g = sp.sum(flow)/(BCs[0]-BCs[1])    # macro conductance
#        print('flow: ', sp.sum(flow))
#        print('dBCs: ', BCs[0]-BCs[1])
#        print('K: ', K)
        return g

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

    def calc_rel_permeability(self, **kwargs):
        try:
            g = self._conductance
        except:
            raise Exception('The conceptual effective conductance is not available')

        if g[-1] - g[0] > 0:
            # Wetting percolation case
            K_r = g/g[0]
        else:
            # Drying percolation case
            K_r = g/g[-1]

        self._rel_permeability = K_r
        return K_r

    def calc_corrected_permeability(self, corr, **kwargs):
        self._corr_permeability = self._rel_permeability*corr
        return self._corr_permeability

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
            pores1[~sp.in1d(p1,P)] = p2[~sp.in1d(p1,P)] #pores1 = BC pores
            pores2[~sp.in1d(p1,P)] = p1[~sp.in1d(p1,P)] #pores2 = inner pores
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

    def run2(self, modify=True, **kwargs):
        r"""

        """
        logger.info("Setup " + self.__class__.__name__)
        self.setup(**kwargs)
        if modify:
            self._modify_system(**kwargs)
        self.solve(**kwargs)

#    @profile
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

        # Excluding the non-spanning/surface pores, see also self.solve()
        try:
            surfpores = self._net.pores('deadend')  # surface pores
        except:
            surfpores = self._net.pores('surface')  # surface pores
        outpores = sp.r_[bcpores, surfpores]
        outpores = sp.unique(outpores)  # pores to ignore
        bcpores = outpores
        inpores = sp.setdiff1d(inpores, surfpores)

        self._bcpores = bcpores
#==============================================================================
#         # Exclude deadend-part of spanning pores
#         # This works, but gives bad result!
#         span = self._net.pores('span')
#         mask = sp.in1d(span, inpores)   # exclude bcpores
#         span = span[mask]
#         num_z = self._net.num_neighbors(pores=span)
#         pdead = num_z == 1
#         self._deadspan = span[pdead] # used in solve()
#         bcpores = sp.unique(sp.r_[bcpores, span[pdead]])
#         inpores = sp.setdiff1d(inpores, span[pdead])
#==============================================================================

        self._inpores = inpores
        A = self.A
        b = self.b

        # Modify rhs b
        bb = A*b    # = 1*Pc_bc (1 in A, and Pc_bc in b, for bcpores)
        self.b = -bb[inpores]
        # Cut BC rows & cols of A
        self.A = self.reshape_mat(A, bcpores)
        if row_scaling:
            self._row_scaling()

    def reshape_mat(self, Amat, id_to_drop):
        '''reshaping matrixwise'''
        C = Amat.tocoo()
        # delete rows & columns
        keep = ~sp.in1d(C.col, id_to_drop)
        C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
        C.col -= id_to_drop.searchsorted(C.col)    # decrement column indices
        C.row -= id_to_drop.searchsorted(C.row)    # decrement column indices
        C._shape = (C.shape[0] - len(id_to_drop), C.shape[1] - len(id_to_drop))
        return C.tocsr()

    def _row_scaling(self, **kwargs):
        r'''
        Modify matrix A by its diagonal, and rhs as well.
        '''
        v = 1./self.A.diagonal()
        nv = len(v)
        self.A = sp.sparse.spdiags(v, 0, nv, nv) * self.A

        b = self.b.T * v
        self.b = sp.reshape(b, self.b.shape)

    def _check_symmetric(self, **kwargs):
        sym = sp.all(self.A-self.A.T==0)
        return sp.all(sym)

    def _check_positive_definite(self, **kwargs):
        pd = sp.all(sp.sparse.linalg.eigsh(self.A, k=100,
                                           return_eigenvectors=0)>0)
        return pd

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
            bc_prop = 'pore.' + self._phase.name +'_bcval_Dirichlet'
            bc_pc = self[bc_prop]
            Dir_pores_vals = bc_pc[self._bcpores]
            self[self._quantity] = 0.0
            self[self._quantity][self._inpores] = self.X
            self[self._quantity][self._bcpores] = Dir_pores_vals
            # deadend-spanning pores
            try:
                deadspan = self._deadspan
                nextodead = self._net.find_neighbor_pores(pores=deadspan, flatten=False)
                self[self._quantity][deadspan] = self[self._quantity][nextodead]
            except:
                pass
            # replace the nan values with the BC pcs
            # as a consequence of excluding the non-spanning pores
            idnan = self._bcpores[sp.isnan(Dir_pores_vals)]
            if sp.any(idnan):
                # Assign pc1 or pc2 to deadend pores (bcpores[idnan])
                inletc = self._net.pores('inlet_clusters')
                outletc = self._net.pores('outlet_clusters')
                ipc = bc_pc[self._net['pore.inlet']].mean()
                opc = bc_pc[self._net['pore.outlet']].mean()
                self[self._quantity][inletc] = ipc
                self[self._quantity][outletc] = opc
        except:
            self[self._quantity] = self.X[self.Ps]
        self[self._quantity+'_abs'] = -self[self._quantity]
        logger.info('Writing the results to ' + '[\'' + self._quantity +
                    '\'] in the ' + self.name + ' algorithm.')

    def store_result(self, Pc=None, sat=None, sat_surf=None, sat_moist=None,
                     w_sat=1, k=None, g=None, k2 = None, dv=None, dv2=None):
        r"""
        Storing Pc, sat, w and k in the object
        """
        if Pc is not None:
            self._Pc = sp.array(Pc)
            self._log_Pc = sp.log10(-self._Pc)
        if sat is not None:
            self._saturation = sp.array(sat)
            self._moisturecontent = self._saturation*w_sat
        if k is not None:
            self._permeability = sp.array(k)
        if k2 is not None:
            self._permeability2 = sp.array(k2)
        if dv is not None:
            self._permeability_vap = sp.array(dv)
        if dv2 is not None:
            self._permeability_vap2 = sp.array(dv2)
        if g is not None:
            self._conductance = sp.array(g)

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

    def return_rate(self, pores=None, throats=None, case='',
                    debug_values=False, pc=0, save_pcrate=False, **kwargs):
        r"""
        Send rate of simulation to phase.
        """
        if pores is None:
            pores = self.Ps
        if throats is None:
            throats = self.Ts

        phase = self._phase
        phase_quantity = self._quantity.replace(phase.name + '_', '')
        if phase_quantity not in phase.props():
            phase[phase_quantity] = sp.nan
            phase[phase_quantity+'_abs'] = sp.nan
        phase[phase_quantity][pores] = self[self._quantity][pores]
        phase[phase_quantity+'_abs'][pores] = self[self._quantity+'_abs'][pores]
        conn_arr = self._net.find_connected_pores(self.Ts)
        dPc = sp.squeeze(sp.diff(self[self._quantity][conn_arr], n=1, axis=1))
        g = self['throat.conductance']
        rate = sp.multiply(g, dPc)
        if 'throat.rate' not in phase.props():
            phase['throat.rate_'+case] = sp.nan
        phase['throat.rate_'+case][throats] = rate[throats]
        phase['throat.delta_pressure_'+case] = dPc
        logger.debug('Results of ' + self.name +
                     ' algorithm have been added to ' + phase.name)
        if save_pcrate:
            pcpores = self[self._quantity]
            phase['pore.pc_'+case+str(int(pc))] = pcpores
            phase['pore.dpc0_'+case+str(int(pc))] = sp.absolute(pcpores - pc)
            phase['throat.rate_'+case+str(int(pc))] = rate
            phase['throat.delta_pressure_'+case+str(int(pc))] = dPc
            phase['throat.conduit_conductance_'+case+str(int(pc))] = g

        if debug_values:
            try:
                phase._values
            except:
                phase._values = {'wetting': {'pc': [], 'rate': {}, 'dp': {},
                                'g_conduit': {}, 'ppressure': {},
                                'poccupy':{}, 'toccupy': {}, 'econd': []},
                                'drying_wetting': {'pc': [], 'rate': {}, 'dp': {},
                                'g_conduit': {}, 'ppressure': {},
                                'poccupy':{}, 'toccupy': {}, 'econd': []},
                                'imbibition': {'pc': [], 'rate': {}, 'dp': {},
                                'g_conduit': {}, 'ppressure': {},
                                'poccupy':{}, 'toccupy': {}, 'econd': []},}
            values = phase._values
            water = self._net._phases[0]

            values[case]['pc'].append(pc)
            values[case]['rate'][pc] = rate
            values[case]['dp'][pc] = dPc
            values[case]['g_conduit'][pc] = g
            values[case]['ppressure'][pc] = self[self._quantity] #phase['pore.pressure_'+case]
            values[case]['poccupy'][pc] = water['pore.occupancy_'+case]
            values[case]['toccupy'][pc] = water['throat.occupancy_'+case]
            values[case]['econd'].append(self._econd[0])



if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
