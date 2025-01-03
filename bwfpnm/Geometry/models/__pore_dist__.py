# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 21:24:59 2015

@author: islah

NB: translation from Staf's net90.F90
"""
import scipy as sp
from scipy import interpolate


def psd(alfa, nr, lrM, lr0, wpc=None, pc=None):
    sigma = 0.07258
    ldr = (lrM - lr0)/nr
    ilist = sp.arange(0, nr)
    lr = lr0 + (ilist + 0.5)*ldr
    r = sp.power(10, lr)
    pcr = -2*sigma/r
#    if pc is not None:
#        pcr = pc
    if wpc is None:
        dwdp = dwdpc(-pcr)      # dw(pc)/dpc
    else:
        f = interpolate.interp1d(pc, wpc)
        w = f(pcr)
        dpc = sp.gradient(pcr)
        dwdp = -sp.gradient(w, dpc)  # dw(pc)/dpc

    # pore volume distribution:
    # fv(r) = dw(pc)/dlog10(pc) = dw(pc)/dpc*dpc/dlog10(pc)
    # = dw/dpc*pc*ln(10)
    pvd = dwdp*pcr*sp.log(10)
#    pvd = movingaverage(pvd, 5)


#    # calculation of integral
#    r11 = sp.zeros_like(r)
#    in_x = sp.zeros_like(r)
#    r11 = sp.power(r, -(alfa + 2))
#    in_x[0] = r11[0] * pvd[0] * ldr
#    for i in range(1, nr):
#        in_x[i] = in_x[i-1] + r11[i]*pvd[i]*ldr
#
#    # cumulative pore radius distribution
#    cum_prd = sp.zeros(sp.size(in_x)+1)
#    cum_prd[1:] = in_x/in_x[-1]
#    lr_prd = sp.r_[lr0, lr]                 # log10(r)

    lr_prd, cum_prd = cum_prd_from_pvd(lr, pvd, alfa, lr0, ldr)
#    print('sum of abs diff lr_prd: ', sum(abs(lr_prd-lr_prd2)))
#    print('sum of abs diff cum_prd: ', sum(abs(cum_prd-cum_prd2)))
#    print('diff cum_prd: ', cum_prd-cum_prd2)

    # normalized cumulative pore volume distribution
    cum_pvd = norm_cum_func(lr, pvd)

    # return (log(r), cum_prd, pvd, cum_pvd)
    return (lr_prd, cum_prd, lr, pvd, cum_pvd)


def cum_prd_from_pvd(lr, pvd, alfa, lr0, ldr):
    r''' Normalized Cumulative pore radius distribution from pvd'''
    nr = sp.size(lr)
    r = sp.power(10, lr)

    # calculation of integral
    r11 = sp.zeros_like(r)
    in_x = sp.zeros_like(r)
    r11 = sp.power(r, -(alfa + 2))
    in_x[0] = r11[0] * pvd[0] * ldr
    for i in range(1, nr):
        in_x[i] = in_x[i-1] + r11[i]*pvd[i]*ldr

    # cumulative pore radius distribution
    cum_prd = sp.zeros(sp.size(in_x)+1)
    cum_prd[1:] = in_x/in_x[-1]
    lr_prd = sp.r_[lr0, lr]                 # log10(r)
    return (lr_prd, cum_prd)


def norm_cum_func(lr, pvd):
    r'''Normalized cumulative pore volume distribution'''
    cum_pvd = sp.integrate.cumtrapz(pvd, x=lr)
    cum_pvd = sp.r_[0, cum_pvd]
    cum_pvd = cum_pvd/max(cum_pvd)
    return cum_pvd


def func_from_cum(x, cum_func):
    dx = sp.gradient(x)
    f = sp.gradient(cum_func, dx)
    return f


def dwdpc(cp):
    reta, retn, retw = sp.zeros(5), sp.zeros(5), sp.zeros(5)
    #	cementpasta
    wsat = 150.
    reta[0] = 1.3e-4
    reta[1] = 2.9e-4
    retn[0] = 1.85
    retn[1] = 3.50
    retw[0] = 0.99
    retw[1] = 0.01
    #	oolieten
    wsat = 310.
    reta[0] = 8.0e-7
    reta[1] = 9.0e-6
    retn[0] = 4.27
    retn[1] = 1.8
    retw[0] = 0.35
    retw[1] = 0.65

    #	calcium silicate
    wsat = 312.09
    reta[0] = 3.4e-7
    reta[1] = 2.0e-6
    reta[2] = 8.0e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 2.00
    retw[0] = 0.21
    retw[1] = 0.50
    retw[2] = 0.056
    #baksteen
    wsat = 159.72
    reta[0] = 27.6326e-6
    reta[1] = 17.1529e-6
    reta[2] = 383.5042e-6
    retn[0] = 1.6914
    retn[1] = 4.4567
    retn[2] = 1.3614
    retw[0] = 0.4544
    retw[1] = 0.4913
    retw[2] = 0.0542
    #baksteen
    wsat = 240.
    reta[0] = 9.9e-6
    reta[1] = 2.4e-5
    reta[2] = 3.835e-4
    retn[0] = 1.638
    retn[1] = 4.0
    retn[2] = 1.5
    retw[0] = 0.165
    retw[1] = 0.7808
    retw[2] = 0.0542
    #alcium silicate
    wsat = 238.992
    reta[0] = 3.4e-7
    reta[1] = 2.5e-6
    reta[2] = 4.50e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 4.00
    retw[0] = 0.274151
    retw[1] = 0.652741
    retw[2] = 0.073107
    #alcium silicate
    wsat = 312.0
    reta[0] = 3.4e-7
    reta[1] = 2.5e-6
    reta[2] = 8.0e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 4.00
    retw[0] = 0.21
    retw[1] = 0.50
    # ictief
    wsat = 100.0
    reta[0] = 100.0e-6
    retn[0] = 1.8
    retw[0] = 1.0
    # ellular concrete
    wsat = 800.0
    reta[0] = 4.0e-9
    reta[1] = 9.2e-7
    reta[2] = 5.0e-6
    reta[3] = 6.0e-5
    reta[4] = 3.0e-3
    retn[0] = 2.0
    retn[1] = 2.3
    retn[2] = 4.0
    retn[3] = 3.5
    retn[4] = 3.5
    retw[0] = 0.013
    retw[1] = 0.25
    retw[2] = 0.04
    retw[3] = 0.07
    retw[4] = 0.627
    # 	calcium silicate plate
    wsat = 894.0
    reta[0] = 2.46e-6
    reta[1] = 2.0e-6
    reta[2] = 2.47e-6
    retn[0] = 1.5
    retn[1] = 2.0
    retn[2] = 4.80
    retw[0] = 0.05
    retw[1] = 0.1
    retw[2] = 0.758

    # this is the derivative of w(pc) w.r.t. pc. w(pc) based on van Genuchten
    # Question: m_i = 1-1/ni, but here is -1+1/ni
    dum = 0.0
    for i in range(sp.size(reta)):
        dum4 = sp.power(cp*reta[i], retn[i])
        dum3 = 1.0 + dum4
        dum2 = sp.log(dum3)
        dum1 = -1.0 + 1.0/retn[i]
        dum0 = sp.exp(dum2*dum1)
        dum = dum + retw[i]*dum0*dum1*dum4*retn[i]/(dum3*cp)

    return dum*wsat


def wpc(cp):
    reta, retn, retw = sp.zeros(5), sp.zeros(5), sp.zeros(5)
    #	cementpasta
    wsat = 150.
    reta[0] = 1.3e-4
    reta[1] = 2.9e-4
    retn[0] = 1.85
    retn[1] = 3.50
    retw[0] = 0.99
    retw[1] = 0.01
    #	oolieten
    wsat = 310.
    reta[0] = 8.0e-7
    reta[1] = 9.0e-6
    retn[0] = 4.27
    retn[1] = 1.8
    retw[0] = 0.35
    retw[1] = 0.65

    #	calcium silicate
    wsat = 312.09
    reta[0] = 3.4e-7
    reta[1] = 2.0e-6
    reta[2] = 8.0e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 2.00
    retw[0] = 0.21
    retw[1] = 0.50
    retw[2] = 0.056
    #baksteen
    wsat = 159.72
    reta[0] = 27.6326e-6
    reta[1] = 17.1529e-6
    reta[2] = 383.5042e-6
    retn[0] = 1.6914
    retn[1] = 4.4567
    retn[2] = 1.3614
    retw[0] = 0.4544
    retw[1] = 0.4913
    retw[2] = 0.0542
    #baksteen
    wsat = 240.
    reta[0] = 9.9e-6
    reta[1] = 2.4e-5
    reta[2] = 3.835e-4
    retn[0] = 1.638
    retn[1] = 4.0
    retn[2] = 1.5
    retw[0] = 0.165
    retw[1] = 0.7808
    retw[2] = 0.0542
    #alcium silicate
    wsat = 238.992
    reta[0] = 3.4e-7
    reta[1] = 2.5e-6
    reta[2] = 4.50e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 4.00
    retw[0] = 0.274151
    retw[1] = 0.652741
    retw[2] = 0.073107
    #alcium silicate
    wsat = 312.0
    reta[0] = 3.4e-7
    reta[1] = 2.5e-6
    reta[2] = 8.0e-5
    retn[0] = 1.85
    retn[1] = 1.60
    retn[2] = 4.00
    retw[0] = 0.21
    retw[1] = 0.50
    # ictief
    wsat = 100.0
    reta[0] = 100.0e-6
    retn[0] = 1.8
    retw[0] = 1.0
    # ellular concrete
    wsat = 800.0
    reta[0] = 4.0e-9
    reta[1] = 9.2e-7
    reta[2] = 5.0e-6
    reta[3] = 6.0e-5
    reta[4] = 3.0e-3
    retn[0] = 2.0
    retn[1] = 2.3
    retn[2] = 4.0
    retn[3] = 3.5
    retn[4] = 3.5
    retw[0] = 0.013
    retw[1] = 0.25
    retw[2] = 0.04
    retw[3] = 0.07
    retw[4] = 0.627
    # 	calcium silicate plate
    wsat = 894.0
    reta[0] = 2.46e-6
    reta[1] = 2.0e-6
    reta[2] = 2.47e-6
    retn[0] = 1.5
    retn[1] = 2.0
    retn[2] = 4.80
    retw[0] = 0.05
    retw[1] = 0.1
    retw[2] = 0.758

    # this is the derivative of w(pc) w.r.t. pc. w(pc) based on van Genuchten
    # Question: m_i = 1-1/ni, but here is -1+1/ni
    dum = 0.0
    for i in range(sp.size(reta)):
        dum4 = sp.power(cp*reta[i], retn[i])
        dum3 = 1.0 + dum4
        dum2 = sp.log(dum3)
        dum1 = -1.0 + 1.0/retn[i]
        dum0 = sp.exp(dum2*dum1)
        dum = dum + retw[i]*dum0

    return dum*wsat


def psd_rand(nelem, nr, cum_prd, lr, alfa):
#    lr = rvs(nelem, lr, cum_prd)
#    r = sp.power(10, lr)

    r = sp.zeros(nelem)
    for ielem in range(0, nelem):
        random = sp.rand()
        # generate a random radius
        r[ielem] = _radius(cum_prd, lr, random)
    ii = sp.sum(r < sp.power(10, -4.5))
    print('[psd_rand()] Fraction of pores with radii < 1e-4.5: ', ii/nelem)

#    # permeabilities
#    rho_l = 1e3
#    mu = 1e-3
#    dl = rho_l/mu/8
#
#    temp = 273.16 + 25.0
#    Rv = 462.0
#    dfact = 1.0e3 * Rv * temp
#    delta_air = 8.69e-10*sp.power(temp, 1.81)
#    delta = delta_air*_psat(temp)/(dfact*temp*Rv)

    vol, lt = sp.zeros(nelem), sp.zeros(nelem)
#    perm0, permv = sp.zeros(nelem), sp.zeros(nelem)
#    voltot = 0.0
    for ielem in range(0, nelem):
        r1 = r[ielem]
        lt[ielem] = sp.power(r1, alfa)   # L(r)=r^alfa
#        a1 = sp.pi*r1*r1             # A = pi.r^2
#        vol[ielem] = a1*lt[ielem]
#        voltot = voltot + vol[ielem]
#        # g_l = (rho_l.r^2/8mu)(A/L) with dl=rho_l/mu/8
#        perm0[ielem] = dl*r1*r1*a1/lt[ielem]
#        # g_v = k_v.(A/L)
#        permv[ielem] = delta*a1/lt[ielem]
    perm0, permv = 0, 0
    return (r, lt, perm0, permv)


def _radius(fr, lr, fx):
    r'''
    Arguments
    ---------
        fx:   random number
    '''
    nr = sp.size(fr) + 1
    for ir in range(0, nr):
        df1 = fr[ir] - fx
        df2 = fr[ir+1] - fx
        if (df1*df2 < 0.0):
            df = fr[ir+1] - fr[ir]
            dr = lr[ir+1] - lr[ir]
            lrc = lr[ir] - df1*dr/df
            return sp.power(10, lrc)


def rvs(size, x, freq):
    # --- normal distributed random numbers ---
    rand_numbers = sp.random.uniform(low=0.0, high=1.0, size=size)
    # --- cdf ---
#    cdf = sp.cumsum(freq)/sp.sum(freq)
    cdf = freq
    # --- inverse cdf ---
    # find X s.t. P(X<=x) = y for each y in rand_numbers
    xr = sp.zeros(size)
    for i, y in enumerate(rand_numbers):
        ind = len(cdf[cdf <= y]) - 1
        dx = (y - cdf[ind])/(cdf[ind+1] - cdf[ind]) * (x[ind+1] - x[ind])
        xr[i] = x[ind] + dx
    return xr


def add_micropores(x1, vmicropor, cum_prd, lr, alfa):
    x = lr
    x = sp.array(x)
    xmin = min(x)
    Lx = xmin - x1
    xx = x - Lx
    xxmin = xx.min()
    Lxx = sp.absolute(xxmin-xx.max())
    xx = xxmin + (xx - xxmin)*Lx/Lxx

    r, lt = [], []
    vmicro = 0
    while vmicro < vmicropor:
        random = sp.rand()
        r.append(_radius(cum_prd, xx, random))
        lt.append(sp.power(r[-1], alfa))
        vmicro += sp.pi*r[-1]**2*lt[-1]
    return (r, lt)


def _psat(t):
    psat = sp.exp(65.8094 - 7066.27/t - 5.976*sp.log(t))
    return psat


def extent_psd(x, y, xmin, xmax, nx=100):
    from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev
#    spl = InterpolatedUnivariateSpline(x, y, k=1)
    tck = splrep(x,y)
    xnew = sp.random.uniform(xmin, xmax, nx)
    xnew.sort()
#    ynew = spl(xnew)
    ynew = splev(xnew, tck)
    return (xnew, ynew)


def extent_cumpsd(y0, y1, y2, y, x, x1, x2, rate, c, weightx=None, weighty=None):
    r'''To rescale the x or y of cumulative pore volume/radius distribution
    x in log10 scale
    '''
    y1 = 0
    # x1 <= x <= x0
    m = (y2-y0)/(y1-y0)
    yy = m*(y-y0) + y0

    # x2 <= x <= x1
    dx = sp.mean(sp.diff(x))
    xx_ext = sp.arange(x2, x1, dx)
    xmean = sp.average([x1, x2])
    yy_ext = gen_sigmoid(xx_ext[1:], xmean, y1, y2, c, rate, 1)

    xx = sp.r_[x2, xx_ext[1:], x]
    yy = sp.r_[0, yy_ext, yy]
    return (xx, yy)


def extent_refpvd(x1, x, y, micropor=0.2, macropor=1):
    x = sp.array(x)
    xmin = min(x)
    Lx = xmin - x1
    xx = x - Lx
    xxmin = xx.min()
    Lxx = sp.absolute(xxmin-xx.max())
    xx = xxmin + (xx - xxmin)*Lx/Lxx

    weighty = micropor/(macropor)
    yy = y*weighty

    xx = sp.r_[xx, x[1:]]
    yy = sp.r_[yy, y[1:]]
    yy = movingaverage(yy, 10)
    return (xx, yy)


def gen_sigmoid(x, xmean, a, b, c, d, e):
    r''' Generalized sigmoid function
    xmean: x s.t. y(x) = (b-a)/2
    a, b:      minimum, maximum asymptote
    d:         growth rate / curve steepness, >> 1 discontinue
    c:         <1 shift curve to left, >1 to right
    e:         steepness of the left part of the curve'''
    y = a + (b-a)/(1+c*sp.exp(-d*(x-xmean)))**1/e
    return y


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = sp.average(values, weights=weights)
    variance = sp.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, sp.sqrt(variance))


def movingaverage(interval, window_size):
    r''' Moving average to smoothen data'''
    window = sp.ones(int(window_size))/float(window_size)
    return sp.convolve(interval, window, 'same')


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    alfa = -1
    nr = 200
    lr1 = -3
    lr0 = -6
    lr_prd, cum_prd, lr, pvd, cum_pvd = psd(alfa, nr, lr1, lr0)
#    plt.figure('fv')
#    plt.plot(lr, pvd, 'o')
#    plt.figure('fr')
#    plt.plot(lr_prd, cum_prd, 'o')

    # 1: extend cumulative function
    x1 = min(lr_prd)
    x2 = sp.log10(1e-9)
    y0, y1 = max(cum_prd), min(cum_prd)
    y2 = 5/(20+5)
    rate = 3
    c = 0.5
    lr_prd1, cum_prd1 = extent_cumpsd(y0, y1, y2, cum_prd, lr_prd, x1, x2, rate, c)
    frxx = movingaverage(cum_prd1, 10)
    cum_prd1[:-10] = frxx[:-10]

    x1 = min(lr)
    x2 = sp.log10(1e-9)
    y0, y1 = max(cum_pvd), min(cum_pvd)
    y2 = 5/(20+5)
    rate = 3
    c = 0.5
    lr1, cum_pvd1 = extent_cumpsd(y0, y1, y2, cum_pvd, lr, x1, x2, rate, c)#,
#                               weightx=[1,1,1,1], weighty=[1,1,1,1])
    frx3 = movingaverage(cum_pvd1, 10)
    cum_pvd1[:-10] = frx3[:-10]

    plt.figure('cum prd1')
    plt.plot(lr_prd, cum_prd, 'o', lr_prd1, cum_prd1, '^')#, lr_prd1, frxx, 's')
    plt.figure('cum pvd1')
    plt.plot(lr, cum_pvd, 'o', lr1, cum_pvd1, '^')#, lr1, frx3, 's')

    pvd1 = func_from_cum(lr1, cum_pvd1)
    prd1 = func_from_cum(lr_prd1, cum_prd1)
    prd = func_from_cum(lr_prd, cum_prd)

    plt.figure('pvd1')
    plt.plot(lr, pvd/max(pvd), 'o', lr1, pvd1, '^')
    plt.figure('prd1')
    plt.plot(lr_prd, prd/max(prd), 'o', lr_prd1, prd1, '^')

    # 2nd: extend with scaled reference pvd
    lr2, pvd2 = extent_refpvd(x2, lr, pvd, micropor=0.2, macropor=1)

    cum_pvd2 = norm_cum_func(lr2, pvd2)

    nr2 = sp.size(lr2)
    ldr = (max(lr2)-min(lr2))/nr2
    r2 = sp.power(10, lr2)
    lr_prd2, cum_prd2 = cum_prd_from_pvd(lr2, pvd2, alfa, min(lr2), ldr)
    prd2 = func_from_cum(lr_prd2[1:], cum_prd2[1:])

    plt.figure('pvd2')
    plt.plot(lr, pvd, 'o', lr2, pvd2, '^')
    plt.figure('Prd2')
    plt.plot(lr_prd, prd, 'o', lr_prd2[1:], prd2, '^')
    plt.figure('cum pvd2')
    plt.plot(lr, cum_pvd, 'o', lr2, cum_pvd2, '^')
    plt.figure('Cum prd2')
    plt.plot(lr_prd, cum_prd, 'o', lr_prd2, cum_prd2, '^')

    # adding micropores to macropores s.t. vmicropor/vmat = micropor
    nelem = 1000
    rmacro, lmacro, perm0, permv = psd_rand(nelem, nr, cum_prd, lr_prd, alfa)
    vmacro = sp.pi*rmacro**2*lmacro
    vmacro = vmacro.sum()

    micropor = 2/100
    macropor = 20/100
    vmicro = vmacro*micropor/macropor
    rmicro, lmicro = add_micropores(x2, vmicro, cum_prd, lr_prd, alfa)

    rgabung = sp.r_[rmicro, rmacro]
    plt.figure('Micropores added to macropores')
    plt.hist(sp.log10(rgabung), 50, histtype='bar', normed=1,  lw=0)




#    nelem = 400
#    r, lt, perm0, permv = psd_rand(nelem, nr, cum_prd, lr, alfa)
#    plt.figure('hist')
#    plt.hist(r, 100)

#    lrr = sp.log10(r)
#    plt.figure('perm')
#    plt.plot(lrr, perm0, 'ro', lrr, permv, 'bs')

#    ldr = (lr1 - lr0)/nr
#    ilist = sp.arange(0, nr)
#    lr = lr0 + (ilist + 0.5)*ldr
#    rr = sp.power(10, lr)
#    sigma = 0.07258
#    pcr = -2*sigma/rr
#    lpcr = sp.log10(-pcr)
#
##    lpc = sp.arange(1, 9.5, 0.5)
#    pc = pcr
#    lpc = sp.log10(-pc)
#    w = wpc(-pc)
#
#    f = interpolate.interp1d(pc, w)
#    wr = f(pcr)
#    plt.figure('w(pc)')
#    plt.plot(lpc, w, 'r-', lpcr, wr, 'bs')
#
#    dwdpr = dwdpc(-pcr)
#    dpc = sp.gradient(pc)
#    dwdp = sp.gradient(w, dpc)
#    plt.figure('dw/dpc')
#    plt.plot(lpc, dwdp, 'r-', lpcr, -dwdpr, 'bs')
#
#    dlpc = sp.gradient(lpc)
#    pv = sp.gradient(w, dlpc)
#    plt.figure('dw/dlogpc')
#    plt.plot(lpc, pv, 'r-', lpcr, pvd, 'bs')
#
#    lr2, fr2, pv2 = psd(alfa, nr, lr1, lr0, w, pc)
#    plt.figure('fr')
#    plt.plot(lr2, fr2, '-')
