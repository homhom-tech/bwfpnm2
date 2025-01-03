# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 14:36:18 2014

@author: islah

Methods: Pore Volume Distribution, f_V(r) = dw/dlog(r) = - dw/dlog(Pc(r))
pore_vol_dist(w,x, xtype='log_pc')
pore_rad_dist(f_V_logr, r, dlength = 'unity')
pore_length(r, L_M, alpha=-1)
van_genuchten(w_sat, pc=None, w_res=0, mat='ceramic')

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps

from bwfpnm.Geometry.models import __wpc__ as wpc
from bwfpnm.Geometry.models import __paper_data__ as pd
from bwfpnm.Utilities.math_func import interpolation

__author__ = '\n'.join(['Muhammad Islahuddin '
                        '<islah.islahuddin@bwk.kuleuven.be>'])

__all__ = ['pore_vol_dist',
           'pore_rad_dist',
           'pore_rad_dist_logpc',
           'pore_length',
           'van_genuchten_unimodal',
           'van_genuchten_multimodal',
           'paperbased_radii_2d',
           'random_from_dist']


def rvs(size, x, freq):
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


def pore_vol_dist(w, x, xtype='log_pc'):
    """Pore volume distribution \n
    Just calculating the derivative of w w.r.t. x, dw/dx, with additional
    argument xtype: \n
    f_V(r) = dw/dlog(r) = - dw/dlog(Pc(r)). central difference for internal
    grid points"""
    # --- if xtype === 'log_pc':
    y = [-(w[1]-w[0])/(x[1]-x[0])]
    y.extend([-((w[i+2]-w[i+1])/(x[i+2]-p) + (w[i+1]-w[i])/(p-x[i]))/2
              for i, p in enumerate(x[1:-1])])
    y.append(-(w[-1]-w[-2])/(x[-1]-x[-2]))

    if xtype in ['r_log', 'r']:
        y = -y
    return list(y)


def pore_rad_dist(f_V_logr, r, L=None):
    """Pore radius distribution \n
    calculated from pore volume distribution \n
    f_R(r) = (dw/dlog10(r))/(rho_l.pi.r^2.L)"""
    densityofwater = 1e+3
    if L is None:
        L = 1e-2

    y = f_V_logr/(densityofwater*np.pi*L*r**2)
    return list(y)


def pore_rad_dist_logpc(S, pc_log, porosity_total, r, L=None):
    """Pore radius distribution \n
    calculated from scratch \n
    f_R(r) = porosity*(-dS/dlog10(pc))/(rho_l.pi.r^2.L)"""
    if L is None:
        L = 1e-2

    dSdlogpc = np.array(pore_vol_dist(S, pc_log))
    y = porosity_total*(dSdlogpc)/(np.pi*L*r**2)
    return list(y)


def pore_length(r, L_M, alpha=-1):
    """Calculate length distribution of a given radius distribution
    alpha: pore shape ratio"""
    r_M = np.max(r)
    L = L_M*(r/r_M)**alpha
    return L


def van_genuchten_unimodal(w_sat, pc=None, w_res=0, mat='ceramic'):
    """Calculate water retention curve w(pc) for certain material,
    mat=['ceramic','silicate']"""
    if mat == 'ceramic':
        a = 2e-5
        n = 2.5
        m = 1-1/n
    elif mat == 'silicate':
        a = 2e-6
        n = 1.5
        m = 1-1/n

    if pc is None:
        pc = wpc.generate_pc()

    effective_saturation = 1/(1+(a*np.abs(pc))**n)**m
    y = w_res + (w_sat - w_res)*effective_saturation
    return y


def van_genuchten_multimodal(mat, pc=None, w_res=0):
    """Calculate water retention curve w(pc) for certain material,
    mat=['ceramic','silicate']"""

    if pc is None:
        pc = wpc.generate_pc()

    effective_saturation = 0*pc
    n = mat.exponent_n
    c = mat.parameter_c
    gama = mat.weight_factor
    for i, ga in enumerate(gama):
        m = -(1-1/n[i])
        effective_saturation += ga*(1+(c[i]*np.abs(pc))**n[i])**m

    try:
        if mat.case == 'Wetting':
            w_sat = mat.w_cap
        elif mat.case == 'Drainage':
            w_sat = mat.w_sat
    except:
        pass
    # --- temporary assumption ---
    w_sat = mat.w_sat
    w = w_res + (w_sat - w_res)*effective_saturation
    return w


def random_from_dist(x, y, n_pore=100, alpha=-1):
    """ --- Generate random data from a given distribution ---\n
    return: random data(r, L)"""
    y = y[::-1]     # reverse the position order
    x = x[::-1]
    a, b = np.min(x), np.max(x)
    yy = y/sum(y)   # normalization

    random_rlog = rvs(n_pore, x, y)

#    custm = Cont_Dist(a, b, x, yy)
#    random_rlog = custm.rvs(size=n_pore)
    random_rlog[random_rlog < a] = a
    random_rlog[random_rlog > b] = b

    random_r = np.power(10., random_rlog)
    random_r.sort()
    random_rlog.sort()
    random_L = pore_length(random_r, L_M=np.max(random_r)*10, alpha=alpha)

    random_data = dict(r=random_r, L=random_L, r_log=random_rlog)
    return random_data


def plot_poredata(r_log, r_lognew, rj_log, pc_log,
                  f_V_logpc, f_V_logpcnew, w, f_V_r, f_V_logr_new,
                  f_R_logr, f_R_logrj_normal, f_R_logpc, material, case):
    """ ================= Plots ======================== """
    plt.figure('Plot Poredata')       # Pore volume distribution
    plt.clf()
    plt.subplot(221)
    plt.plot(r_log, f_V_logpc, 'o', r_lognew, f_V_logpcnew)
    plt.title(r'Pore Volume Distribution $f_V(r)$')
    plt.xlabel('log(r) [m]')
    plt.ylabel(r'$f_V(r)$')

    plt.subplot(222)
    plt.plot(pc_log, w, pc_log, f_V_logpc)
    plt.title('{} in {}'.format(case, material))
    plt.xlabel('log(pc) [Pa]')
    plt.ylabel(r'w & $f_V(r)$')

    plt.subplot(223)
#    plt.plot(r_log,f_R_logr,r_log,f_R_logpc)
    plt.plot(r_log, f_R_logpc)
    plt.title(r'Pore Radius Distribution $f_R(r)$')
    plt.xlabel('log(r) [m]')
    plt.ylabel(r'$f_R$')
#    plt.plot(r_log, f_V_r)
#    plt.xlabel('log(r) [m]')
#    plt.ylabel(r'$dw/dr$')

    plt.subplot(224)
    plt.title(r'Normalized Pore Radius Distribution $f_j(r)$')
    plt.plot(rj_log, f_R_logrj_normal)
    plt.xlabel('log(r) [m]')
    plt.ylabel(r'$f_j([r_j,r_M])$')
    plt.show()
#    plt.plot(r_log, f_V_logr_new, r_log,f_V_logpcnew)
#    plt.xlabel('log(r) [m]')
#    plt.ylabel(r'$dw/dlog_{10}(r)$')


def paperbased_radii_2d(material='Ceramicbrick', case='adsorption',
                        n_pore=10, fig_plot=False):
    """Generate random pore sizes based on known material properties:
    ceramic and silicate bricks"""

    pc = -np.logspace(3, 9, n_pore)
    pc_log = np.log10(-pc)
    r = wpc.pc_from_r(pc)
    r_log = np.log10(r)

    if material == 'Ceramicbrick':
        mat = pd.ceramicbrick_carmeliet99()
#        mat = paper_data.ceramicbrick_carmeliet01(case)
    elif material == 'Calciumsilicate':
        mat = pd.calciumsilicate_carmeliet99()
#        mat = pd.calciumsilicate_carmeliet01(case)

    w = van_genuchten_multimodal(mat, pc)

    f_V = np.array(pore_vol_dist(w, np.log10(-pc)))
    L = pore_length(r, L_M=np.max(r)*6, alpha=mat.pore_shape_ratio)
    f_R = np.array(pore_rad_dist(f_V, r, L))

    poredistribution = dict(porosity_total=mat.open_porosity, w=w, pc=pc,
                            pc_log=pc_log, r_log=r_log, f_R=f_R, f_V=f_V)

    # --- calculate random data based on random radius from f_R_normal ---
    rand_data = random_from_dist(r_log, f_R, n_pore, mat.pore_shape_ratio)
#    rand_data = random_from_dist(r_log, f_R, n_pore, 1.5)
    rand_pore_vol = np.pi*rand_data['L']*rand_data['r']**2
    rand_pore_frac = rand_pore_vol/rand_pore_vol.sum()

    rand_data.update({'pore_frac': rand_pore_frac})

    if fig_plot:
        w_sat = mat.w_sat
        S = w/w_sat
        porosity_total = mat.open_porosity
        r_min = mat.radius_min
        r_max = mat.radius_max
        pore_vol = np.pi*L*r**2
        pore_frac = pore_vol/pore_vol.sum()

        rand_pclog = np.log10(-wpc.pc_from_r(rand_data['r']))
        rand_f_R = interpolation(r_log, f_R, xnew=rand_data['r_log'])
        rand_f_V = interpolation(r_log, f_V, xnew=rand_data['r_log'])
        rand_w = interpolation(pc_log, w, xnew=rand_pclog)

        plt.rc('text',usetex=True)
        font = {'family':'serif','size':16}
        plt.rc('font',**font)
        plt.rc('legend',**{'fontsize':14})
        plt.rc('axes',**{'labelweight':'normal', 'labelsize':14})

        f, axarr = plt.subplots(2, 2)
        rand_rlog = rand_data['r_log']
        axarr[0, 0].plot(r_log, f_V,rand_rlog, rand_f_V, 'ro')
        axarr[0, 0].legend(['given data','random data'], loc='upper left')
        axarr[0, 0].set_title('Pore Volume Distribution',fontsize=14)
        axarr[0, 0].set_xlabel('$log(r)$')

        axarr[0, 1].plot(pc_log, w, pc_log, f_V, rand_pclog,rand_w,'ro')
        axarr[0, 1].set_title('Moisture Content vs Pore Volume Distribution',fontsize=14)
        axarr[0, 1].set_xlabel('$log_{10}(-P_c)$')

        axarr[1, 0].plot(r_log, f_R, rand_rlog, rand_f_R,'ro')
        axarr[1, 0].set_title('Pore Radius Distribution',fontsize=14)
        axarr[1, 0].set_xlabel('$log(r)$')

        #%% --- Different scale of observation: Normalized f_R ---
        # ----- not used yet: for later multiscale approach
        rj_list=  np.array([r_min*10**i for i in range(5) if r_min*10**i < r_max])
        r_log_min = np.log10(rj_list)      # list of rj_log
        smntr = zip(r_log,f_R)

        for i,rj in enumerate(rj_list):
            rj_log = np.array([x for x,y in smntr if x > r_log_min[i]])
            f_Rj = np.array([y for x,y in smntr if x > r_log_min[i]])

            porosity_j = porosity_total*(1-S[r == rj])
            nj_pores = simps(f_Rj,-rj_log)
            f_Rj_normal = f_Rj/nj_pores   # integral(f_Rj_normal,-rj_log)=1
    #            axarr[1, 1].plot(rj_log, f_Rj_normal)
    #        axarr[1, 1].set_title('Normalized Pore Radius Distribution')
    #        axarr[1, 1].set_xlabel('log(r)')

        axarr[1, 1].plot(r_log,pore_frac[::-1].cumsum()[::-1],
                        rand_rlog,rand_pore_frac.cumsum(),'ro')
        axarr[1, 1].set_title('Pore Fraction',fontsize=14)
        axarr[1, 1].set_xlabel('$log(r)$')

        plt.suptitle('{} in {}'.format(case, material), fontsize=16)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        plt.pause(1)

    return rand_data, poredistribution


if __name__ == '__main__':
    mat = {0: 'Ceramicbrick', 1: 'Calciumsilicate'}
    cases = {0: 'Wetting', 1: 'Drainage'}
#    rand_data, poredistribution = generate_paperbased_radii(material = mat[0],
#                                                            case = cases[0],
#                                                            n_pore = 100,
#                                                            fig_save = None)
    ni, nj = (20, 20)
    n_throat = (ni-1)*nj + (nj-1)*ni
    rand_data, poredistribution = paperbased_radii_2d(mat[1], cases[0],
                                                      n_throat)
    plt.figure()
    plt.plot(poredistribution['pc_log'], poredistribution['w'])
    plt.xlabel('log(-$p_c$)')
    plt.figure()
    plt.plot(poredistribution['r_log'], poredistribution['f_V'])
    plt.xlabel('log(r)')
    plt.figure()
    plt.plot(poredistribution['r_log'], poredistribution['f_R'])
    plt.xlabel('log(r)')
    plt.figure()
    plt.hist(rand_data['r'])
