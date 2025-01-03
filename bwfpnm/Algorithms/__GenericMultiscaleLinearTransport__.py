# -*- coding: utf-8 -*-
"""
===============================================================================
module __GenericLinearTransport__: Class for solving linear transport processes
for multiscale cases
===============================================================================
"""
import scipy as sp
import scipy.sparse.linalg as sprslin
from OpenPNM.Algorithms import GenericLinearTransport
from OpenPNM.Base import logging
try:
    from pyamg import smoothed_aggregation_solver, rootnode_solver,\
                        ruge_stuben_solver, solve
    from pyamg.aggregation import adaptive_sa_solver
except:
    pass
import warnings
import time
from functools import wraps
try:
    from petsc4py import PETSc
    nopetsc = False
except:
    nopetsc = True
from scipy.sparse.linalg import eigsh
logger = logging.getLogger(__name__)

try:
    solvers_dict = {'rs': ruge_stuben_solver, 'classic': ruge_stuben_solver,
                    'sa': smoothed_aggregation_solver,
                    'rootnode': rootnode_solver, 'ra': rootnode_solver,
                    'blackbox': solve, 'asa': adaptive_sa_solver}
except:
    solvers_dict = {}


class GenericMultiscaleLinearTransport(GenericLinearTransport):
    r"""
    Customized GenericLinearTransport with Algebraic Multigrid solvers.
    alg.run(amg)

    Available options for amg solvers:
    'ruge_stuben_solver': a classic amg, RS
    'smoothed_aggregation_solver': SA
    'rootnode_solver': combination of RS and SA solvers
    'blackbox': for arbitrary A

    """
    def __init__(self, phase=None, **kwargs):
        super().__init__(phase, **kwargs)

    def _fn_timer(function):
        @wraps(function)
        def function_timer(*args, **kwargs):
            t0 = time.clock()
            result = function(*args, **kwargs)
            t1 = time.clock()
            print ("Total time running %s: %s seconds" %
                   (function.func_name, str(t1-t0))
                   )
            return result
        return function_timer

    def _find_backbone(self, Mat, threshold=1e-4):
        r'''Identify linearly dependent rows from the matrix.
        Such rows correspond to the dead-end part of the spanning clusters!

        Use this function when need it.
        https://stackoverflow.com/questions/28816627/how-to-find-linearly-independent-rows-from-a-matrix

        Return indices of the linearly independent columns

        Created: 15 Mar 2018
        Status: problematic in determining the threshold
        '''
        M = Mat.copy()
        pivot = sp.absolute(M.data).min()
        M.data /= pivot

        norm = sprslin.norm
        absolute = sp.absolute

        dim = sp.shape(Mat)[0]
        value = sp.sparse.lil_matrix((dim, dim))
        index = []
        index.append(0) #without loss of generality we pick the first column as linearly dependent
        MMt = M.dot(M.T)
        Mnorm = [norm(M[:,i]) for i in range(dim)]
        Mnorm = sp.array(Mnorm)
        j = 0 #therefore the second index is simply 0
        for i in range(dim): #loop over the columns
            if i != j: #if the two columns are not the same
                # CS: dot(u,v) <= norm(u)*norm(v)
                # u, v linearly dependen/parallel <=> dot(u,v) = norm(u)*norm(v)
                value[i,j] = absolute(MMt[i,j] - Mnorm[i] * Mnorm[j])
                if value[i,j] > threshold:  # record the lin. indenpendent instead
                    index.append(i) # index of lin. ind. pore is saved
                    j = i #j is refreshed
        index = sp.array(index)
        return index

    def _do_one_inner_iteration(self, A, b, amg=None, x0=None,
                                maxiter=5000, tol=1e-14,
                                strength=('symmetric', {'theta': 0.03}),
                                CF='RS', agg='standard',
                                smooth=('gauss-seidel', {'omega': 4.0/3.0}),
                                max_coarse=1000, max_levels=10,
                                accel=None, cycle='F', res=None,
                                save_grids=False, keep=False,
                                parallel=False, umfpack=True, econd=False,
                                backbone=False, **kwargs):
        r"""
        This method solves AX = b and returns the result to the corresponding
        algorithm.

        Solver:
        -------
        solvers: a list of str solvers
            'direct' --> The default solver is SuperLU
                (included in the scipy distribution), which can solve real or
                complex linear systems in both single and double precisions.
                It is automatically replaced by UMFPACK, ifavailable.
            'cg', 'gmres', etc...
            'pyamg_rs', 'pyamg_sa', 'pyamg_ra': pyamg
            'petsc_gamg', 'petsc_ml', 'petsc_boomeramg': PETSc

        direct solver
            amg=None (default), and iterative_solver=None (default)

        iterative solver
            amg=None, and iterative_solver = 'cg', 'gmres'

        amg iterative solver
            amg='rs', 'classic', 'sa', 'rootnode', 'blackbox'
        """

        logger.info('Solving sparse AX = b')

        if A is None:
            A = self.A
        if b is None:
            b = self.b

        if res is None:
            res = []
        try:
            x0 = sp.ravel(sp.ones_like(b))*sp.ravel(x0)
        except:
            x0 = None

#        print('Matrix A, min = {}, max = {}, average = {}'.format(A.data.min(),
#              A.data.max(), A.data.mean()))
        #%% Estimate the condition number
        # -------------------------------
        if econd:
            evals_large, evecs_large = eigsh(A, 1, which='LM')
            evals_small, evecs_small = eigsh(A, 1, sigma=0, which='LM')
            econd = evals_large/evals_small
            self._econd = econd
        #%% Identify the backbone part of the spanning cluster, & exclude the deadend part
        if backbone:
            boneids = self._find_backbone(A)
        #%% Solving phase
        if parallel and not nopetsc:
            X = self.parallel_solver(A, b, x0, **kwargs)
            self._solver = 'PETSc'
            return X
        elif parallel and nopetsc:
            raise Exception('PETSc is not available')

        if amg is not None:
#            try:
#                # skip construction phase if perm_err is True
##                print('pc {}'.format(x0[0]))
##                raise(RuntimeError)
##                print('skipping the raise')
#                X = self._amg_ml.solve(b, x0=x0, tol=tol, accel=accel,
#                                       residuals=res, maxiter=maxiter, cycle=cycle)
#
#            except:
            try:
                self._solver = 'pyamg-' + amg + '-' + accel
            except:
                self._solver = 'pyamg-' + amg
            try:

                amg = solvers_dict[amg.lower()]
            except:
                logger.error('No such amg method is available! Set amg=None.')

            logger.info('Iterative solver: AMG, %s'.format(amg.__name__))
            if sp.any(save_grids):
                keep = True

            if amg.__name__ == 'ruge_stuben_solver':
                ml = amg(A,
                         strength=strength,
                         CF=CF,
                         presmoother=('gauss_seidel', {'sweep': 'symmetric'}),
                         postsmoother=('gauss_seidel', {'sweep': 'symmetric'}),
                         max_levels=max_levels, max_coarse=max_coarse, keep=keep)
            elif amg.__name__ == 'smoothed_aggregation_solver':
                ml = amg(A, B=None, BH=None, symmetry='hermitian',
                         strength=strength,
                         aggregate=agg,
                         smooth=smooth,
                         presmoother=('block_gauss_seidel',
                                     {'sweep': 'symmetric'}),
                         postsmoother=('block_gauss_seidel',
                                      {'sweep': 'symmetric'}),
                         improve_candidates=[('block_gauss_seidel',
                                            {'sweep': 'symmetric',
                                            'iterations': 4}), None],
                         max_levels=max_levels, max_coarse = max_coarse,
                         diagonal_dominance=False, keep=keep)
            elif amg.__name__ == 'rootnode_solver':
                ml = amg(A, B=None, BH=None, symmetry='hermitian',
                        strength=strength,
                        aggregate=agg,
                        smooth=smooth,
                        presmoother=('block_gauss_seidel',
                                     {'sweep': 'symmetric'}),
                        postsmoother=('block_gauss_seidel',
                                      {'sweep': 'symmetric'}),
                        improve_candidates=[('block_gauss_seidel',
                                            {'sweep': 'symmetric',
                                             'iterations': 4}), None],
                        max_levels=max_levels, max_coarse = max_coarse,
                        diagonal_dominance=False, keep=keep)

            try:
                X = ml.solve(b, x0=x0, tol=tol, accel=accel,
                             residuals=res, maxiter=maxiter, cycle=cycle)
            except:
                X = solve(A, b, x0=x0, verb=True, tol=tol)

            if sp.size(res) > maxiter and maxiter > 1:
                message = 'The AMG iteration does not converge for the given\
                            tolerance within maxiter'
                logger.warning(message)
                warnings.warn(message)

        elif self._iterative_solver is None:
#            cA = sp.absolute(self.A.data).min()
#            cb = sp.absolute(self.b).max()
#            cA, cb = 1, 1
            try:

                X = sprslin.spsolve(A, b, use_umfpack=umfpack)#(A/cA, b/cb)

            except:
                X = sp.zeros_like(self.b)
#            X *= cb/cA
            self._solver = 'direct'
            if umfpack:
                self._solver += '_umfpack'
            else:
                self._solver += '_superlu'

        else:
            if self._iterative_solver not in ['cg', 'gmres']:
                raise Exception('GenericLinearTransport does not support the' +
                                ' requested iterative solver!')
            self._solver = 'iterative-' + self._iterative_solver
            params = kwargs.copy()
            solver_params = ['x0', 'tol', 'maxiter', 'xtype', 'M', 'callback']
            [params.pop(item, None) for item in kwargs.keys()
             if item not in solver_params]
            tol = kwargs.get('tol')
            if tol is None:
                tol = 1e-20
            params['tol'] = tol
            if self._iterative_solver == 'cg':
                result = sprslin.cg(A, b, **params)
            elif self._iterative_solver == 'gmres':
                result = sprslin.gmres(A, b, **params)
            X = result[0]
            self._iterative_solver_info = result[1]

        try:
            if sp.ceil(x0[0]) in sp.ceil(save_grids):
                try:
                    self._amg_ml
                    ml
                except:
                    try:
                        self._amg_ml = ml
                    except:
                        pass

                if sp.any(save_grids):
                    self._create_multigrid_array(sp.ceil(x0[0]), X, **kwargs)
        except:
            pass

        return X

    def _create_multigrid_array(self, pc, X, **kwargs):
        r'''
        Save to the physic instance the multigrid structure as well as
        the corresponding conduit conductance and occupancy data
        '''
        for phase in self._net._phases:
            if phase.name == 'water':
                case = self.name.split('_')[-1]
                water = phase
        try:
            inpores = self._inpores         # if modify_mat = True
        except:
            inpores = self._net.pores()     # if otherwise
        self._net['pore.internal'] = self.tomask(inpores)

        pc = case + '_' + str(abs(int(pc)))
        prop_level = 'pore.levels_' + pc
        prop_cond = 'throat.conduit conductance_' + pc
        prop_occ = 'occupancy_' + pc
        prop_sol = 'pore.pc_' + self._solver + '_' + pc
        obj_ins = self._phase._physics[0]
        obj_ins[prop_level] = 0

        obj_ins[prop_cond] = self['throat.conductance']
        obj_ins['pore.'+prop_occ] = water['pore.occupancy_'+case]
        obj_ins['throat.'+prop_occ] = water['throat.occupancy_'+case]
        # Save the Pc field
        # Removing the additional super pore variables from the results
        try:
            Dir_pores_vals = self['pore.' + self._phase.name +
                                          '_bcval_Dirichlet'][self._bcpores]
            obj_ins[prop_sol] = 0.0
            obj_ins[prop_sol][self._inpores] = X
            obj_ins[prop_sol][self._bcpores] = Dir_pores_vals
        except:
            obj_ins[prop_sol] = X[self.Ps]

        try:
            index0 = sp.where(self._amg_ml.levels[0].splitting == 1)[0]
            obj_ins[prop_level][inpores[index0]] = 1
            splittings = [index0]
            for i, level in enumerate(self._amg_ml.levels[1:-1]):
                index = sp.where(level.splitting == 1)[0]
                splittings.append(index)
                for coarse in splittings[-2::-1]:
                    index = coarse[index]
                obj_ins[prop_level][inpores[index]] = i+2

        except:     # no multigrid found
            message = 'No multigrid found, only 1 level is available'
            logger.warning(message)
            warnings.warn(message)

    def parallel_solver(self, A, b, x, **kwargs):
        r'''
        Solve in parallel using PETSc package
        '''
        Ap = self._create_pmat(A)
        bp = self._create_pvec(b)
        xp = self._create_pvec(x)
        ksp = self._create_psolver(Ap, bp, xp, **kwargs)
        return ksp.getSolution().getArray()

    def _create_pmat(self, M):
        r'''
        Create a sparse csr  PETSc matrix from a csr matrix M
        '''
        A = PETSc.Mat()
        A.create(PETSc.COMM_WORLD)
        A.createAIJ(size=M.shape, csr=(M.indptr, M.indices, M.data))
        A.setType('mpiaij')
        A.assemblyBegin()
        A.assemblyEnd()
        return A

    def _create_pvec(self, array):
        r'''
        Create a sparse csr  PETSc matrix
        '''
        return PETSc.Vec().createWithArray(array)

    def _create_psolver(self, A, b, x=None, solve=True, **kwargs):
        r'''
        Create linear solver

        isolver: KSP linear solver objects
            'preonly', 'gmres', 'cg', 'richardson', 'chebyshev',
            'bicg', 'fgmres', 'dgmres', 'gcr', 'bcgs', 'cgs', 'tfqmr',
            'tcqmr', 'cr', 'lsqr'
        prec: PC preconditioners
            'none', 'gamg', 'lu', 'ilu', 'hypre', 'ml', 'jacobi',
            'bjacobi', 'sor', 'eisenstat', 'icc', 'asm', 'gasm',
            'bddc', 'ksp', 'composite', 'cholesky', 'shell'

            External solvers (KSPType='preonly', PCType='lu', ):
            'matlab', 'mumps', 'superlu', 'superlu_dist', 'umfpack',
            'essl', 'lusol'


        To check:
        - solver: ksp.getType()
        - preconditioner: ksp.getPC().getType()
        '''
        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        # set options
        opts = PETSc.Options()
        PETSc.Options.Set('ksp_monitor')
        opts.Set('ksp_vecmonitor')
        opts_dict = {}

        # set the preconditioner and its options (PC)
        self._petsc_pc_options(kwargs['pc_type'], opts_dict, **kwargs)
        for key, val in opts_dict.items():
            opts[key] = val
        # set linear operator & solve
        ksp.setOperators(A)
        ksp.view()
        ksp.setFromOptions()
        ksp.view()
        if solve:
            ksp.solve(b, x)
        return ksp

    def _petsc_pc_options(self, prec, opts_dict, **kwargs):
        opts_dict.update(self._ksp_defaults)
        opts_dict.update(self._defamg[prec])
        for key in kwargs.keys():
            if key in opts_dict.keys():
                opts_dict[key] = kwargs[key]

    def _compare_np_pet_mats(self, N, P):
        acsr = P.getValuesCSR()
        assert sp.all(acsr[2]==N.data)
        assert sp.all(acsr[1]==N.indices)
        assert sp.all(acsr[0]==N.indptr)
        print('Both matrices are identic!')

    def _compare_np_pet_vecs(self, n, p):
        assert sp.all(n.flatten()==p.getArray().flatten())
        print('Both vectors are identic!')

    # ref: http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCML.html
#            ksp.pc.setHYPREType('boomeramg')
#        or: mpiexec -n 2 ./ex2 -pc_type hypre -pc_hypre_type boomeramg -help |grep pc_hypre_boomeramg_
    _ksp_defaults = {'ksp_type': 'gmres',
                     'ksp_gmres_restart': 20,
                     'ksp_initial_guess_nonzero': True,
                     'ksp_max_it': 5000,
                     'ksp_rtol': 1e-10,
                     'pc_type': 'gamg',
                     'pc_mg_cycles': 2,
                     'pc_mg_smoothup': 1,
                     'pc_mg_smoothdown': 1,
                     'pc_mg_levels': 10,
                     'log_view': 1,
                     'log_summary': 1,
                     'ksp_view':1}
    _gamg_defaults = {'pc_gamg_type': 'classical',
                      'pc_gamg_smooths': 1,
                      'pc_gamg_agg_nsmooths': 1, # or 0: unsmoothed aggregation
                      'pc_gamg_threshold': 0.02,
                      'pc_gamg_repartition': False,
                      'pc_gamg_reuse_interpolation': False,
                      'pc_gamg_asm_use_agg': False,
                      'pc_gamg_use_parallel_coarse_grid_solver': True,
                      'pc_gamg_process_eq_limit': 100,
                      'pc_gamg_coarse_eq_limit': 100,
                      'pc_gamg_sym_graph': True,
                      'pc_gamg_square_graph': 1,
                      'mg_levels_ksp_type': 'chebyshev',
                      'mg_levels_pc_type': 'sor',
                      'mg_levels_ksp_max_it': 2,
                      'mg_coarse_ksp_type': 'preonly',
                      'mg_coarse_pc_type': 'lu'
                      }
    _hypre_defaults = {'pc_hypre_type': 'boomeramg',
                      'pc_hypre_boomeramg_cycle_type': 'w',
                      'pc_hypre_boomeramg_max_levels': 25,
                      'pc_hypre_boomeramg_max_iter': 1,
                      'pc_hypre_boomeramg_rtol': 0.0,
                      'pc_hypre_boomeramg_tol': 0,
                      'pc_hypre_boomeramg_truncfactor': 0,
                      'pc_hypre_boomeramg_P_max': 0,
                      'pc_hypre_boomeramg_agg_nl': 10,
                      'pc_hypre_boomeramg_agg_num_paths': 1,
                      'pc_hypre_boomeramg_strong_threshold': 0.25,
                      'pc_hypre_boomeramg_max_row_sum': 0.9,
                      'pc_hypre_boomeramg_nodal_coarsen': 0,
                      'pc_hypre_boomeramg_vec_interp_variant': 0,
                      'pc_hypre_boomeramg_grid_sweeps_all': 1,
                      'pc_hypre_boomeramg_grid_sweeps_down': 1,
                      'pc_hypre_boomeramg_grid_sweeps_up': 1,
                      'pc_hypre_boomeramg_grid_sweeps_coarse': 1,
                      'pc_hypre_boomeramg_smooth_type': 'Schwarz-smoothers',
                      'pc_hypre_boomeramg_smooth_num_levels': 25,
                      'pc_hypre_boomeramg_eu_level': 0,
                      'pc_hypre_boomeramg_eu_droptolerance': 0.,
                      'pc_hypre_boomeramg_eu_bj': 0,
                      'pc_hypre_boomeramg_relax_type_all': 'symmetric-SOR/Jacobi',
                      'pc_hypre_boomeramg_relax_type_down': 'symmetric-SOR/Jacobi',
                      'pc_hypre_boomeramg_relax_type_up': 'symmetric-SOR/Jacobi',
                      'pc_hypre_boomeramg_relax_type_coarse': 'Gaussian-elimination',
                      'pc_hypre_boomeramg_no_CF': 0,
                      'pc_hypre_boomeramg_coarsen_type': 'Falgout',
                      'pc_hypre_boomeramg_interp_type': 'classical',
                      'pc_hypre_boomeramg_relax_weight_all': 1,
                      'pc_hypre_boomeramg_relax_weight_level': (1,1),
                      'pc_hypre_boomeramg_outer_relax_weight_all': 1,
                      'pc_hypre_boomeramg_outer_relax_weight_level': (1,1),
                      'pc_hypre_boomeramg_measure_type':'local',
                      'pc_hypre_boomeramg_nodal_relaxation': 0,
                      'pc_hypre_boomeramg_print_statistics': 1,
                      'pc_hypre_boomeramg_print_debug': 1}
    _ml_defaults = {'pc_ml_maxNlevels': 10,
                      'pc_ml_maxCoarseSize': 500,
                      'pc_ml_Threshold': 0,
                      'pc_ml_CoarsenScheme': 'Uncoupled',
                      'pc_ml_PrintLevel': 0,
                      'pc_ml_DampingFactor': 1.33333,
                      'pc_ml_SpectralNormScheme_Anorm': False,
                      'pc_ml_Symmetrize': True,
                      'pc_ml_BlockScaling': False,
                      'pc_ml_nullspace': 'AUTO',
                      'pc_ml_EnergyMinimization': 0,
                      'pc_ml_reuse_interpolation': False,
                      'pc_ml_KeepAggInfo': False,
                      'pc_ml_Reusable': True,
                      'pc_ml_OldHierarchy': False,
                      'pc_ml_repartition': False}
    _defamg ={'gamg': _gamg_defaults,
              'hypre': _hypre_defaults,
              'ml': _ml_defaults}
