# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 11:08:54 2014

@author: islah
"""
#import os
#import sys
#sys.path.append(os.path.abspath(os.path.dirname(__file__) + '/' + '../..'))

import scipy.stats as ss
import numpy as np
import scipy.sparse.linalg as spla
from numpy import linalg as lan
from scipy import interpolate
#from pysparse.sparse import spmatrix
#from pysparse.sparse.pysparseMatrix import PysparseMatrix
#from pysparse.direct.pysparseSuperLU import PysparseSuperLUSolver
#from pysparse.direct.pysparseUmfpack import PysparseUmfpackSolver
#from pysparse.precon import precon
#from pysparse.itsolvers import krylov

__author__ = '\n'.join(['Muhammad Islahuddin '
                        '<islah.islahuddin@bwk.kuleuven.be>'])

__all__ = ['Cont_Dist',
           'interpolation',
           'check_matrix_regularity',
           'pysparse_matrix',
           'scipy_solver',
           'pysparse_solver']


class Cont_Dist(ss.rv_continuous):
    """ ---------- FAST!! -------------
    A class for building a specific continuous distribution.

    The distribution is defined by a set of data (x,y).

    Example
    ---------

    """

    def __init__(self, a, b, sb_x, Freqs):
        """
        a   :   left boundary
        b   :   right boundary
        sb_x:   array of x
        Freqs:  f(x)
        """
        super(Cont_Dist, self).__init__(a=a, b=b)
        self.sumbu_x = sb_x
        self.Freqs = Freqs

    def _cdf(self, x):
        return [self.do_cdf(i) for i in x]

    def do_cdf(self, x):
        if x < self.a:
            return 0.0
        if x >= self.b:
            return 1.0

        freq = self.Freqs[self.sumbu_x <= x]
        v = freq.sum()

        ind = len(freq) - 1
        v += (self.Freqs[ind+1] - self.Freqs[ind]) * (x - self.sumbu_x[ind]) /\
             (self.sumbu_x[ind+1] - self.sumbu_x[ind])
        return v


def rvs(size, x, freq):
    r'''
    Generate a 'size' numbers of samples based on a given distribution (x, freq)

    Parameters
    ----------
    size : int
        Number of samples to be drawn

    x, freq : list, list
        supplied distribution
    '''
    # --- normal distributed random numbers ---
    rand_numbers = np.random.uniform(low=0.0, high=1.0, size=size)
    # --- cdf ---
    cdf = np.cumsum(freq)/np.sum(freq)
    # --- inverse cdf ---
    # find X s.t. P(X<=x) = y for each y in rand_numbers
    xr = np.zeros(size)
    for i, y in enumerate(rand_numbers):
        ind = len(cdf[cdf <= y]) - 1
        dx = (y - cdf[ind])/(cdf[ind+1] - cdf[ind]) * (x[ind+1] - x[ind])
        xr[i] = x[ind] + dx
    return xr


def interpolation(x, y, Npoint=100, xnew=None, method='cubicspline'):
    """interpolate from sorted x,y to xnew,ynew using cubic-spline

    Required packages:
    from scipy import interpolate

    arguments:
    ----------
    x, y            : original data
    Npoint=100      : number points of new x to be generated
    xnew=None       : specified new x
    method='cubicspline'  : ['linear','nearest', 'zero', 'slinear',
                             'quadratic', 'cubic']
    return:
    ---------
    xnew, ynew      : if xnew is None
    ynew            : otherwise
    """
    # sorting the x axis
    sortfirst = list(zip(x, y))
    sortfirst.sort()
    x = [n for (n, m) in sortfirst]
    y = [m for (n, m) in sortfirst]
    xmin, xmax = np.min(x), np.max(x)

    # generate new & finer data
    if xnew is None:
        xnew = np.linspace(xmin, xmax, Npoint)

    if method == 'cubicspline':
        tck = interpolate.splrep(x, y, s=0)
        ynew = interpolate.splev(xnew, tck, der=0)
    elif method in {'linear', 'nearest', 'zero', 'slinear',
                    'quadratic', 'cubic'}:
        tck = interpolate.interp1d(x, y, kind=method)
        ynew = tck(xnew)

    if xnew is None:
        return xnew, ynew
    else:
        return ynew


#%% --- MATRICES ---
def check_matrix_regularity(A):
    """Check the regularity of a matrix of a linear system
    Input: A - a sparse (tridiagonal) matrix"""

    A_dense = A.todense()
    eigen_values = lan.eigvals(A_dense)
    determinant = lan.det(A_dense)
    condition_number = lan.cond(A_dense)
    try:
        lan.cholesky(-A.todense())
    except:
        posdef = 'no'
#        print '(-A) is not positive definite'
    else:
        posdef = 'yes'
#        print '(-A) is positive definite'

    try:
        lan.inv(A.todense())
    except:
        invertible = 'no'
#        print 'A is singular or not square'
    else:
        invertible = 'yes'
#        print 'A is invertible'

    regularity = dict(eigen=eigen_values, det=determinant,
                      cond=condition_number, posdef=posdef,
                      invertible=invertible)
    return regularity


def scipy_solver(A, b, solver=('direct', 'superilu'), precond='ilu',
                 x0=None, tol=1e-10):
    """ --------- Scipy's solver for sparse system -------
    Solve A.x = b using method indicated in solver
    """
    A = A.astype(np.float64)
    if solver[0] == 'direct':
        info = solver[0]
        if solver[1] == 'umfpack':
            spla.use_solver(assumeSortedIndices=False)
            sol = spla.spsolve(A, b, use_umfpack=True)
        elif solver[1] == 'superlu':
            sol = spla.spsolve(A, b, use_umfpack=False)
        elif solver[1] == 'superilu':
            try:
                solve = spla.spilu(A.tocsc(), drop_tol=1e-8)
            except:
                solve = spla.spilu(A.tocsc(), drop_tol=1e-4)
                print('warning! RunTimeError when running spilu \
                      with drop_tol = 1e-8')

            sol = solve.solve(b)

    elif solver[0] == 'iterative':
        if precond == 'ilu':
            solve = spla.spilu(A.tocsc(), drop_tol=1e-8)

            def matvec(xx):
                return solve.solve(xx)

            def rmatvec(xx):
                return solve.solve(xx, 'T')

            M = spla.LinearOperator(A.shape,
                                    matvec=matvec, rmatvec=rmatvec)
#            x0 = solve.solve(b)
        else:
            M = None

        if solver[1] == 'bicg':
            sol, info = spla.bicg(A, b, x0=x0, tol=tol, M=M)
        elif solver[1] == 'cgs':
            sol, info = spla.cgs(A, b, x0=x0, tol=tol, M=M)
        elif solver[1] == 'gmres':
            sol, info = spla.gmres(A, b, x0=x0, tol=tol, M=M)
        elif solver[1] == 'lgmres':
            sol, info = spla.lgmres(A, b, x0=x0, tol=tol, M=M)
        elif solver[1] == 'qmr':
            sol, info = spla.qmr(A, b, x0=x0, tol=tol, M1=M, M2=M)

    relres = np.linalg.norm(A*sol-b)/np.linalg.norm(b)

    return sol, relres, info


def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a-b) <= np.maximum(rel_tol * np.maximum(abs(a), abs(b)), abs_tol)




if __name__ == '__main__':
    x = np.arange(10)
    y = max(x) - x
    size = 1e+4
    xr = rvs(size, x, y)
    import matplotlib.pyplot as plt
    count, bins, ignored = plt.hist(xr, 15, normed=True)
    plt.plot(bins, np.ones_like(bins), linewidth=2, color='r')
    plt.show()
