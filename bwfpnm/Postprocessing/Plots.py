import scipy as _sp
import matplotlib.pylab as _plt
from bwfpnm.Utilities.IO import save_data_csv


def setting_fig():
    xmin, xmax = _plt.xlim()
    ymin, ymax = _plt.ylim()
    _plt.xlim(xmin, xmax)
    _plt.ylim(ymin, ymax)
    _plt.axvline(xmin, color='k', linewidth=3)
    _plt.axhline(ymin, color='k', linewidth=3)
    _plt.axvline(xmax, color='k', linewidth=3)
    _plt.axhline(ymax, color='k', linewidth=3)
    _plt.xticks(fontsize=15)
    _plt.yticks(fontsize=15)

def profiles(network, values=None, bins=[10,10,10]):
    r'''
    Compute the profiles for the property of interest and plots it in all
    three dimensions

    Parameters
    ----------
    network : OpenPNM Network object

    values : array_like, optional
        The pore property values to be plotted as a profile

    bins : int or list of ints, optional
        The number of bins to divide the domain into for averaging.

    Notes
    -----
    Either propname or values can be sent, but not both

    '''
    fig = _plt.figure()
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    ax = [ax1,ax2,ax3]
    xlab = ['x_coordinate','y_coordinate','z_coordinate']
    for n in [0,1,2]:
        n_min, n_max = [_sp.amin(network['pore.coords'][:,n]), _sp.amax(network['pore.coords'][:,n])]
        steps = _sp.linspace(n_min,n_max,bins[n]+1,endpoint=True)
        vals = _sp.zeros_like(steps)
        for i in range(0,len(steps)-1):
            temp = (network['pore.coords'][:,n] > steps[i])*(network['pore.coords'][:,n] <= steps[i+1])
            vals[i] = _sp.mean(values[temp])
        yaxis = vals[:-1]
        xaxis = (steps[:-1] + (steps[1]-steps[0])/2)/n_max
        ax[n].plot(xaxis,yaxis,'bo-')
        ax[n].set_xlabel(xlab[n])
        ax[n].set_ylabel('Slice Value')


def distributions(net,
                  fig = None,
                  throat_diameter='throat.diameter',
                  pore_diameter='pore.diameter',
                  throat_length='throat.length',
                  exclude_boundaries=True,
                  geom_list=None,
                  logscale=True,
                  histtype='list',
                  normed=False):
    r"""
    Plot a montage of key network size distribution histograms

    Parameters
    ----------
    net : OpenPNM Network Object
    The network for which the graphs are desired

    """
    if fig is None:
        fig = _plt.figure()

    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.4)

    if geom_list is not None:
        include_pores = [False]*net.num_pores()
        include_throats = [False]*net.num_throats()
        for geom in geom_list:
            include_pores = include_pores | net["pore."+geom]
            include_throats = include_throats | net["throat."+geom]
    else:
        include_pores = net["pore.all"]
        include_throats = net["throat.all"]
    pores = net.pores()[include_pores]
    throats = net.throats()[include_throats]

    Dpore = _sp.log10(net[pore_diameter][pores])
    Dthroat = _sp.log10(net[throat_diameter][throats])
    Lthroat = _sp.log10(net[throat_length][throats])
    Zcoord = net.num_neighbors(pores,flatten=False)

    if histtype=='list':
        ax1 = fig.add_subplot(221)
        ax1.hist(Dpore,25,facecolor='green')
        ax1.set_xlabel('Pore Diameter [log(m)]', fontdict={'fontsize':18})
        ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #    ax1.locator_params(tight=True, nbins=3)
        max_xticks = 3
        xloc = _plt.MaxNLocator(max_xticks)
        ax1.xaxis.set_major_locator(xloc)
        setting_fig()

        ax2 = fig.add_subplot(222)
    #    x = net.num_neighbors(pores,flatten=False)
        bins = _sp.arange(0, Zcoord.max()+1)
        ax2.hist(Zcoord, bins, facecolor='yellow')
        ax2.set_xlabel('Coordination Number', fontdict={'fontsize':18})
        ax2.set_ylabel('Frequency', fontdict={'fontsize':18})
        xloc = _plt.MaxNLocator(max_xticks)
        ax2.xaxis.set_major_locator(xloc)
        setting_fig()

        ax3 = fig.add_subplot(223)
        ax3.hist(Dthroat,25,facecolor='blue')
        ax3.set_xlabel('Throat Diameter [log(m)]', fontdict={'fontsize':18})
        ax3.set_ylabel('Frequency', fontdict={'fontsize':18})
        ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        xloc = _plt.MaxNLocator(max_xticks)
        ax3.xaxis.set_major_locator(xloc)
        setting_fig()

        ax4 = fig.add_subplot(224)
        ax4.hist(Lthroat,25,facecolor='red')
        ax4.set_xlabel('Throat Length [log(m)]', fontdict={'fontsize':18})
        ax4.set_ylabel('Frequency', fontdict={'fontsize':18})
        ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        xloc = _plt.MaxNLocator(max_xticks)
        ax4.xaxis.set_major_locator(xloc)
        setting_fig()
        fig.show()

    elif histtype == 'shift':
        ax1 = fig.add_subplot(121)
        ax1.hist((Dpore, Dthroat, Lthroat), 25, normed=normed, edgecolor='none')
        ax1.set_xlabel('Diameter [log(m)]', fontdict={'fontsize':18})
        ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        max_xticks = 3
        xloc = _plt.MaxNLocator(max_xticks)
        ax1.xaxis.set_major_locator(xloc)
        setting_fig()

        ax2 = fig.add_subplot(122)
        bins = _sp.arange(0, Zcoord.max()+1)
        ax2.hist(Zcoord, bins, facecolor='yellow', normed=normed)
        ax2.set_xlabel('Coordination Number', fontdict={'fontsize':18})
        ax2.set_ylabel('Frequency', fontdict={'fontsize':18})
        xloc = _plt.MaxNLocator(max_xticks)
        ax2.xaxis.set_major_locator(xloc)
        setting_fig()

    elif histtype=='size':
        ax1 = fig.add_subplot(111)
        xmin = _sp.amin(_sp.r_[Dpore, Dthroat, Lthroat])
        xmax = _sp.amax(_sp.r_[Dpore, Dthroat, Lthroat])
#        bins = range(int(_sp.floor(xmin)), int(_sp.ceil(xmax)))
        bins = _sp.linspace(int(_sp.floor(xmin)), int(_sp.ceil(xmax)), 20)
#        bins = _sp.linspace(-7, -3.5, 20)
        ax1.hist((Dpore, Dthroat, Lthroat), bins, normed=normed, edgecolor='none')
        ax1.set_xlabel('Diameter [log(m)]', fontdict={'fontsize':18})
        ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        max_xticks = 5
        xloc = _plt.MaxNLocator(max_xticks)
        ax1.xaxis.set_major_locator(xloc)
        _plt.yscale('log')
        setting_fig()

    elif histtype=='coord':
        ax2 = fig.add_subplot(111)
        bins = _sp.arange(0, Zcoord.max()+1)
        ax2.hist(Zcoord, bins, facecolor='yellow', normed=normed)
        ax2.set_xlabel('Coordination Number', fontdict={'fontsize':18})
        ax2.set_ylabel('Frequency', fontdict={'fontsize':18})
        xloc = _plt.MaxNLocator(max_xticks)
        ax2.xaxis.set_major_locator(xloc)
        setting_fig()


def psd1graph(net,
             fig = None,
             bins=20,
             cum=False,
             normed=False,
             throat_diameter='throat.diameter',
             pore_diameter='pore.diameter',
             exclude_boundaries=True,
             geom_list=None,
             label=('pore', 'throat'),
             color=('blue', 'red'),
             save_csv=False, filename='psd'):
    r"""
    Plot a montage of key network size distribution histograms

    Parameters
    ----------
    net : OpenPNM Network Object
    The network for which the graphs are desired

    """
    if fig is None:
        fig = _plt.figure(figsize=(15,10))

    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.4)

    if geom_list is not None:
        ipore = [False]*net.num_pores()
        ithroat = [False]*net.num_throats()
        for geom in geom_list:
            include_pores = ipore | net["pore."+geom]
            include_throats = ithroat | net["throat."+geom]

    else:
        include_pores = net["pore.all"]
        include_throats = net["throat.all"]

    pores = net.pores()[include_pores]
    throats = net.throats()[include_throats]

    Dpore = _sp.log10(net[pore_diameter][pores])
    Dthroat = _sp.log10(net[throat_diameter][throats])

    if _sp.isscalar(bins):
        binss = _sp.linspace(_sp.amin([Dpore.min(), Dthroat.min()]),
                             _sp.amax([Dpore.max(), Dthroat.max()]), bins)
    else:
        binss = bins

    ax1 = fig.add_subplot(111)
    hist, bins, patches = ax1.hist((Dpore, Dthroat), bins=binss,
                                   cumulative=cum, normed=normed,
                                   label=label, edgecolor='none', color=color)
    ax1.set_xlabel('Diameter [log(m)]', fontdict={'fontsize':18})
    ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#    ax1.locator_params(tight=True, nbins=3)
    max_xticks = 5
    max_yticks = 5
    xloc = _plt.MaxNLocator(max_xticks)
    ax1.xaxis.set_major_locator(xloc)
    yloc = _plt.MaxNLocator(max_yticks)
    ax1.yaxis.set_major_locator(yloc)
    _plt.legend()
    setting_fig()
    fig.show()
    if save_csv:
        data = {'bins': bins, 'pore_diameter':_sp.r_[hist[0],0],
                'throat_diameter':_sp.r_[hist[1],0]}
        save_data_csv(data, filename=filename)
    return (hist, bins)


def psd2nets(net, net2,
             fig = None,
             throat_diameter='throat.diameter',
             pore_diameter='pore.diameter',
             throat_length='throat.length',
             exclude_boundaries=True,
             geom_list=None,
             label=('net1', 'net2'),
             color=('blue', 'red')):
    r"""
    Plot a montage of key network size distribution histograms

    Parameters
    ----------
    net : OpenPNM Network Object
    The network for which the graphs are desired

    """
    if fig is None:
        fig = _plt.figure()

    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.4)

    if geom_list is not None:
        ipore = [False]*net.num_pores()
        ithroat = [False]*net.num_throats()

        for geom in geom_list:
            include_pores = ipore | net["pore."+geom]
            include_throats = ithroat | net["throat."+geom]

        ipore = [False]*net2.num_pores()
        ithroat = [False]*net2.num_throats()
        for geom in geom_list:
            include_pores2 = ipore | net2["pore."+geom]
            include_throats2 = ithroat | net2["throat."+geom]

    else:
        include_pores = net["pore.all"]
        include_throats = net["throat.all"]

        include_pores2 = net2["pore.all"]
        include_throats2 = net2["throat.all"]

    pores = net.pores()[include_pores]
    throats = net.throats()[include_throats]

    pores2 = net2.pores()[include_pores2]
    throats2 = net2.throats()[include_throats2]

    Dpore = _sp.log10(net[pore_diameter][pores])
    Dthroat = _sp.log10(net[throat_diameter][throats])
    Lthroat = _sp.log10(net[throat_length][throats])
    Zcoord = net.num_neighbors(pores,flatten=False)

    Dpore2 = _sp.log10(net2[pore_diameter][pores2])
    Dthroat2 = _sp.log10(net2[throat_diameter][throats2])
    Lthroat2 = _sp.log10(net2[throat_length][throats2])
    Zcoord2 = net2.num_neighbors(pores2,flatten=False)

    binss=_sp.linspace(-6, -3.5, 20)

    ax1 = fig.add_subplot(221)
    ax1.hist((Dpore, Dpore2), binss, label=label,
              edgecolor='none', color=color)
    ax1.set_xlabel('Pore Diameter [log(m)]', fontdict={'fontsize':18})
    ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#    ax1.locator_params(tight=True, nbins=3)
    max_xticks = 5
    max_yticks = 5
    xloc = _plt.MaxNLocator(max_xticks)
    ax1.xaxis.set_major_locator(xloc)
    yloc = _plt.MaxNLocator(max_yticks)
    ax1.yaxis.set_major_locator(yloc)
    _plt.legend()
    setting_fig()

    ax2 = fig.add_subplot(222)
#    x = net.num_neighbors(pores,flatten=False)
    bins = _sp.arange(0, _sp.amax([Zcoord.max(), Zcoord2.max()]))
    ax2.hist((Zcoord, Zcoord2), bins, edgecolor='none', color=color)
    ax2.set_xlabel('Coordination Number', fontdict={'fontsize':18})
    ax2.set_ylabel('Frequency', fontdict={'fontsize':18})
#    ax2.xaxis.set_major_locator(xloc)
    ax2.yaxis.set_major_locator(yloc)
    setting_fig()

    ax3 = fig.add_subplot(223)
    ax3.hist((Dthroat, Dthroat2),binss, edgecolor='none', color=color)
    ax3.set_xlabel('Throat Diameter [log(m)]', fontdict={'fontsize':18})
    ax3.set_ylabel('Frequency', fontdict={'fontsize':18})
    ax3.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax3.xaxis.set_major_locator(xloc)
    ax3.yaxis.set_major_locator(yloc)
    setting_fig()

    ax4 = fig.add_subplot(224)
    ax4.hist((Lthroat, Lthroat2),binss, edgecolor='none', color=color)
    ax4.set_xlabel('Throat Length [log(m)]', fontdict={'fontsize':18})
    ax4.set_ylabel('Frequency', fontdict={'fontsize':18})
    ax4.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax4.xaxis.set_major_locator(xloc)
    ax4.yaxis.set_major_locator(yloc)
    setting_fig()
    fig.show()


def psd2nets2(net, net2,
             fig = None,
             throat_diameter='throat.diameter',
             pore_diameter='pore.diameter',
             throat_length='throat.length',
             exclude_boundaries=True,
             geom_list=None,
             label=('net1', 'net2'),
             color=('blue', 'red')):
    r"""
    Plot a montage of key network size distribution histograms

    Parameters
    ----------
    net : OpenPNM Network Object
    The network for which the graphs are desired

    """
    if fig is None:
        fig = _plt.figure()

    fig.subplots_adjust(hspace = 0.4)
    fig.subplots_adjust(wspace = 0.4)

    if geom_list is not None:
        n = _sp.size(geom_list)
        ipore = [False]*net.num_pores()
        ithroat = [False]*net.num_throats()

        Dpore, Dthroat, Lthroat, Zcoord = [], [], [], []
        include_pores, pores = {}, {}
        include_throats, throats = {}, {}
        for geom in geom_list:
            include_pores[geom] = ipore | net["pore."+geom]
            include_throats[geom] = ithroat | net["throat."+geom]
            pores[geom] = net.pores()[include_pores[geom]]
            throats[geom] = net.throats()[include_throats[geom]]

            Dpore.append(_sp.log10(net[pore_diameter][pores[geom]]))
            Dthroat.append(_sp.log10(net[throat_diameter][throats[geom]]))
            Lthroat.append(_sp.log10(net[throat_length][throats[geom]]))
            Zcoord.append(net.num_neighbors(pores[geom],flatten=False))

        ipore = [False]*net2.num_pores()
        ithroat = [False]*net2.num_throats()
        Dpore2, Dthroat2, Lthroat2, Zcoord2 = [], [], [], []
        include_pores2, pores2 = {}, {}
        include_throats2, throats2 = {}, {}
        for geom in geom_list:
            include_pores2[geom] = ipore | net2["pore."+geom]
            include_throats2[geom] = ithroat | net2["throat."+geom]
            pores2[geom] = net2.pores()[include_pores2[geom]]
            throats2[geom] = net2.throats()[include_throats2[geom]]

            Dpore2.append(_sp.log10(net2[pore_diameter][pores2[geom]]))
            Dthroat2.append(_sp.log10(net2[throat_diameter][throats2[geom]]))
            Lthroat2.append(_sp.log10(net2[throat_length][throats2[geom]]))
            Zcoord2.append(net2.num_neighbors(pores2[geom],flatten=False))


#        include_pores = [False]*net.num_pores()
#        include_throats = [False]*net.num_throats()
#        for geom in geom_list:
#            include_pores = include_pores | net["pore."+geom]
#            include_throats = include_throats | net["throat."+geom]
#
#        include_pores2 = [False]*net2.num_pores()
#        include_throats2 = [False]*net2.num_throats()
#        for geom in geom_list:
#            include_pores2 = include_pores2 | net2["pore."+geom]
#            include_throats2 = include_throats2 | net2["throat."+geom]
    else:
        include_pores = net["pore.all"]
        include_throats = net["throat.all"]

        include_pores2 = net2["pore.all"]
        include_throats2 = net2["throat.all"]

        pores = net.pores()[include_pores]
        throats = net.throats()[include_throats]

        pores2 = net2.pores()[include_pores2]
        throats2 = net2.throats()[include_throats2]

        Dpore = _sp.log10(net[pore_diameter][pores])
        Dthroat = _sp.log10(net[throat_diameter][throats])
        Lthroat = _sp.log10(net[throat_length][throats])
        Zcoord = net.num_neighbors(pores,flatten=False)

        Dpore2 = _sp.log10(net2[pore_diameter][pores2])
        Dthroat2 = _sp.log10(net2[throat_diameter][throats2])
        Lthroat2 = _sp.log10(net2[throat_length][throats2])
        Zcoord2 = net2.num_neighbors(pores2,flatten=False)

    binss=_sp.linspace(-6, -3.5, 22)
    import itertools
    span1 = _sp.r_[Dpore[1], Dpore[2]]
    iso1 = _sp.array(list(itertools.zip_longest(*[span1, Dpore[0]], fillvalue=0)))
    span2 = _sp.r_[Dpore2[1], Dpore2[2]]
    iso2 = _sp.array(list(itertools.zip_longest(*[span2, Dpore2[0]], fillvalue=0)))

    ax1 = fig.add_subplot(111)
    ax1.hist(iso1, binss,
#              label=(label[0], 'iso'),
              edgecolor='none', stacked=True,
              color=('darkblue', 'green'), alpha=0.5)
    ax1.hist(iso2, binss-0.01,
#              label=(label[1], 'iso'),
              edgecolor='none', stacked=True, alpha=0.5,
              color=('darkred', 'yellow'))

    ax1.set_xlabel('Pore Diameter [log(m)]', fontdict={'fontsize':18})
    ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#    ax1.locator_params(tight=True, nbins=3)
    max_xticks = 6
    max_yticks = 5
    xloc = _plt.MaxNLocator(max_xticks)
    ax1.xaxis.set_major_locator(xloc)
    yloc = _plt.MaxNLocator(max_yticks)
    ax1.yaxis.set_major_locator(yloc)
    _plt.legend()
    setting_fig()

#    binss=_sp.linspace(0, 45, 25)
#    import itertools
#    span1 = _sp.r_[Zcoord[1], Zcoord[2]]
#    iso1 = _sp.array(list(itertools.zip_longest(*[span1, Zcoord[0]], fillvalue=0)))
#    span2 = _sp.r_[Zcoord2[1], Zcoord2[2]]
#    iso2 = _sp.array(list(itertools.zip_longest(*[span2, Zcoord2[0]], fillvalue=0)))
#
#    ax1 = fig.add_subplot(111)
#    ax1.hist(iso1, binss,
##              label=(label[0], 'iso'),
#              edgecolor='none', stacked=True,
#              color=('darkblue', 'green'), alpha=0.5)
#    ax1.hist(iso2, binss-0.25,
##              label=(label[1], 'iso'),
#              edgecolor='none', stacked=True, alpha=0.5,
#              color=('darkred', 'yellow'))
#
#    ax1.set_xlabel('Connection Number', fontdict={'fontsize':18})
#    ax1.set_ylabel('Frequency', fontdict={'fontsize':18})
#    ax1.set_xlim([0, 45])
##    ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
##    ax1.locator_params(tight=True, nbins=3)
#    max_xticks = 4
#    max_yticks = 5
#    xloc = _plt.MaxNLocator(max_xticks)
#    ax1.xaxis.set_major_locator(xloc)
#    yloc = _plt.MaxNLocator(max_yticks)
#    ax1.yaxis.set_major_locator(yloc)
#    _plt.legend()
#    setting_fig()


def wetting_curves(inv_alg,
                   Pc='inv_Pc',
                   sat='inv_sat',
                   seq='inv_seq'):
  r"""
  Plot a montage of key saturation plots:
  - Sat(pc)
  - Sat(sim_step)
  - Pc(sat)
  - Pc(sim_step)

  Parameters
  ----------
  inv_alg : OpenPNM Algorithm Object
    The invasion algorithm for which the graphs are desired

  Examples
  --------

  >>> import OpenPNM
  >>> pn = OpenPNM.Network.TestNet()
  >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
  >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
  >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
  >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
  >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=[pn.pores('top')], outlets=pn.pores('bottom'))
       IP algorithm at 0 % completion at 0.0 seconds
       IP algorithm at 20 % completion at 0.0 seconds
       IP algorithm at 40 % completion at 0.0 seconds
       IP algorithm at 60 % completion at 0.0 seconds
       IP algorithm at 100% completion at  0.0  seconds
  >>> OpenPNM.Postprocessing.Plots.drainage_curves(IP,timing='inv_time')

  """
  inv_throats = inv_alg.toindices(inv_alg['throat.'+seq]>0)
  sort_seq = _sp.argsort(inv_alg['throat.'+seq][inv_throats])
  inv_throats = inv_throats[sort_seq]
  tpc_unik = _sp.unique(inv_alg['throat.'+Pc][inv_throats])
  tsat_unik = _sp.unique(inv_alg['throat.'+sat][inv_throats])
  tseq_unik = _sp.unique(inv_alg['throat.'+seq][inv_throats])
  t_unik = _sp.array(list(zip(tpc_unik, tsat_unik, tseq_unik)))

  inv_pores = inv_alg.toindices(inv_alg['pore.'+seq]>0)
  sort_seq = _sp.argsort(inv_alg['pore.'+seq][inv_pores])
  inv_pores = inv_pores[sort_seq]
  ppc_unik = _sp.unique(inv_alg['pore.'+Pc][inv_pores])
  psat_unik = _sp.unique(inv_alg['pore.'+sat][inv_pores])
  pseq_unik = _sp.unique(inv_alg['pore.'+seq][inv_pores])
  p_unik = _sp.array(list(zip(ppc_unik, psat_unik, pseq_unik)))
  unik = _sp.concatenate((t_unik, p_unik))

  unik = unik[unik[:,1].argsort()]

  pt_sat = _sp.unique(unik[:,1])
  pt_seq = _sp.unique(unik[:,2])
  pt_pc = _sp.zeros_like(pt_sat)
  j = 0
  nn = _sp.shape(unik)[0]
  for i in range(0,nn-1):
      if unik[i,2] == unik[i+1,2]:
          pt_pc[j] = max(unik[i,0], unik[i+1,0])
          j += 1
      elif _sp.mod(i, 2)==0:
          pt_pc[j] = unik[i,0]
          j += 1
      elif unik[i,2] != unik[i-1,2]:
          pt_pc[j] = unik[i,0]
          j += 1
  if i+1 == nn-1:
    pt_pc[j] = unik[i+1,0]

#  print 'pc: ', pt_pc
#  print 'sat: ', pt_sat
#  print 'seq: ', pt_seq
#  print 'unik: ', unik

  fig = _plt.figure('Pc Curves', figsize=(13, 10), dpi=80, facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(221)   #left
  ax2 = fig.add_subplot(222)   #middle
  ax4 = fig.add_subplot(223)   #left
  ax5 = fig.add_subplot(224)   #middle

#  ax1.plot(inv_alg['throat.'+Pc][inv_throats],inv_alg['throat.'+sat][inv_throats])
  ax1.plot(pt_pc, pt_sat)
  ax1.set_xlabel('Capillary Pressure (Pa)')
  ax1.set_ylabel('Saturation')
  ax1.set_ylim([0,1])
  ax1.set_xlim([0.99*min(pt_pc), 1.01*max(pt_pc)])
#  ax1.set_xlim([0.99*min(inv_alg['throat.'+Pc][inv_throats]),1.01*max(inv_alg['throat.'+Pc][inv_throats])])

#  ax2.plot(inv_alg['throat.'+seq][inv_throats],inv_alg['throat.'+sat][inv_throats])
  ax2.plot(pt_seq, pt_sat)
  ax2.set_xlabel('Simulation Step')
  ax2.set_ylabel('Saturation')
  ax2.set_ylim([0,1])
  ax2.set_xlim([0,1.01*max(pt_seq)])
#  ax2.set_xlim([0,1.01*max(inv_alg['throat.'+seq][inv_throats])])

#  ax4.plot(inv_alg['throat.'+sat][inv_throats],inv_alg['throat.'+Pc][inv_throats])
  ax4.plot(pt_sat, pt_pc)
  ax4.set_ylabel('Capillary Pressure (Pa)')
  ax4.set_xlabel('Saturation')
  ax4.set_xlim([0,1])
  ax4.set_ylim([0.99*min(pt_pc),1.01*max(pt_pc)])
#  ax4.set_ylim([0.99*min(inv_alg['throat.'+Pc][inv_throats]),1.01*max(inv_alg['throat.'+Pc][inv_throats])])

  ax5.plot(pt_seq,pt_pc)
  ax5.set_xlabel('Simulation Step')
  ax5.set_ylabel('Capillary Pressure (Pa)')
  ax5.set_ylim([0.99*min(pt_pc),1.01*max(pt_pc)])
#  ax5.set_ylim([0.99*min(inv_alg['throat.'+Pc][inv_throats]),1.01*max(inv_alg['throat.'+Pc][inv_throats])])
#  ax5.set_xlim([0,1.01*max(inv_alg['throat.'+seq][inv_throats])])
  ax5.set_xlim([0,1.01*max(pt_seq)])

  fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)
  ax4.grid(True)
  ax5.grid(True)

  fig.show()


def moisture_retention(inv_alg,
                       w_sat,
                       Pc='inv_Pc',
                       sat='inv_sat',
                       seq='inv_seq'):
  r"""
  Plot a montage of key saturation plots:
  - Sat(pc)
  - Sat(sim_step)
  - Pc(sat)
  - Pc(sim_step)

  Parameters
  ----------
  inv_alg : OpenPNM Algorithm Object
    The invasion algorithm for which the graphs are desired

  Examples
  --------

  >>> import OpenPNM
  >>> pn = OpenPNM.Network.TestNet()
  >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
  >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
  >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
  >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
  >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=[pn.pores('top')], outlets=pn.pores('bottom'))
       IP algorithm at 0 % completion at 0.0 seconds
       IP algorithm at 20 % completion at 0.0 seconds
       IP algorithm at 40 % completion at 0.0 seconds
       IP algorithm at 60 % completion at 0.0 seconds
       IP algorithm at 100% completion at  0.0  seconds
  >>> OpenPNM.Postprocessing.Plots.drainage_curves(IP,timing='inv_time')

  """
  # Combining pore & throat data in correct order
  pt_seq, pt_pc, pt_sat = _extract_pt(inv_alg, Pc='inv_Pc', sat='inv_sat',
                                      seq='inv_seq')
  pt_w = pt_sat*w_sat

  fig = _plt.figure('Moisture', figsize=(6, 10), dpi=80, facecolor='w',
                    edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

#  ax1.plot(inv_alg['throat.'+Pc][inv_throats],inv_alg['throat.'+sat][inv_throats])
  ax1.plot(-pt_pc, pt_w)
  ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([0.99*min(-pt_pc), 1.01*max(-pt_pc)])
  ax1.set_xscale('log')

  ax2.plot(pt_pc, pt_sat)
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Saturation [-]')
  ax2.set_ylim([0,1])
  ax2.set_xlim([0.99*min(pt_pc), 1.01*max(pt_pc)])

  fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()


def permeability(sat, pc, eff_permeability):
    r"""
      Plot a montage of key saturation plots:
      - Sat(pc)
      - Sat(sim_step)
      - Pc(sat)
      - Pc(sim_step)

      Parameters
      ----------
      inv_alg : OpenPNM Algorithm Object
        The invasion algorithm for which the graphs are desired

      Examples
      --------

      >>> import OpenPNM
      >>> pn = OpenPNM.Network.TestNet()
      >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
      >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
      >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
      >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
      >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
      >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
      >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=[pn.pores('top')], outlets=pn.pores('bottom'))
           IP algorithm at 0 % completion at 0.0 seconds
           IP algorithm at 20 % completion at 0.0 seconds
           IP algorithm at 40 % completion at 0.0 seconds
           IP algorithm at 60 % completion at 0.0 seconds
           IP algorithm at 100% completion at  0.0  seconds
      >>> OpenPNM.Postprocessing.Plots.drainage_curves(IP,timing='inv_time')

    """
    fig = _plt.figure('Permeability', figsize=(6, 10), dpi=80, facecolor='w',
                      edgecolor='k')
    ax1 = fig.add_subplot(211)   # top
    ax2 = fig.add_subplot(212)   # bottom

#  ax1.plot(inv_alg['throat.'+Pc][inv_throats],inv_alg['throat.'+sat][inv_throats])
    ax1.plot(-pc, eff_permeability)
    ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
    ax1.set_ylabel('Effective Moisture Permeability [kg/msPa]')
    ax1.set_xlim([1.01*min(-pc), 0.99*max(-pc)])
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.plot(sat, eff_permeability)
    ax2.set_xlabel('Saturation [-]')
    ax2.set_ylabel('Effective Moisture Permeability [kg/msPa]')

    fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.1)
    ax1.grid(True)
    ax2.grid(True)

    fig.show()


def moisture_permeability(pc, w, eff_permeability,
                          name='Moisture & Permeability',
                          Pc='inv_Pc',
                          sat='inv_sat',
                          seq='inv_seq',
                          color='k',
                          linestyle='-',
                          marker='o',
                          x1scale='log',
                          y1scale=None,
                          x2scale='log',
                          y2scale='log',
                          legend='adsorption'):
  r"""
  Plot a montage of key saturation plots:
  - Sat(pc)
  - Sat(sim_step)
  - Pc(sat)
  - Pc(sim_step)

  Parameters
  ----------
  inv_alg : OpenPNM Algorithm Object
    The invasion algorithm for which the graphs are desired

  Examples
  --------

  >>> import OpenPNM
  >>> pn = OpenPNM.Network.TestNet()
  >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
  >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
  >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
  >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
  >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=[pn.pores('top')], outlets=pn.pores('bottom'))
       IP algorithm at 0 % completion at 0.0 seconds
       IP algorithm at 20 % completion at 0.0 seconds
       IP algorithm at 40 % completion at 0.0 seconds
       IP algorithm at 60 % completion at 0.0 seconds
       IP algorithm at 100% completion at  0.0  seconds
  >>> OpenPNM.Postprocessing.Plots.drainage_curves(IP,timing='inv_time')

  """
  # Combining pore & throat data in correct order
#  pt_seq, pt_pc, pt_sat = _extract_pt(inv_alg, Pc='inv_Pc', sat='inv_sat',
#                                      seq='inv_seq')
#  pt_w = pt_sat*w_sat

  fig = _plt.figure(name, figsize=(6, 10), dpi=80,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

#  p1, = ax1.plot(-pt_pc, pt_w,
#                 color = color, linestyle = linestyle, marker = marker)

  p1, = ax1.plot(-pc, w,
                 color=color, linestyle=linestyle, marker=marker)
  ax1.set_xlabel('Capillary Pressure [Pa]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([1.01*min(-pc), 0.99*max(-pc)])
  ax1.legend([p1], [legend])
  if x1scale is not None:
      ax1.set_xscale(x1scale)
      ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Log Moisture Content [log(kg/m^3)]')

  p2, = ax2.plot(-pc, eff_permeability,
                 color=color, linestyle=linestyle, marker=marker)
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Effective Moisture Permeability [kg/msPa]')
  ax2.set_xlim([1.01*min(-pc), 0.99*max(-pc)])
  ax2.legend([p2], [legend])
  if x2scale is not None:
      ax2.set_xscale(x2scale)
      ax2.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y2scale is not None:
      ax2.set_yscale(y2scale)
      ax2.set_ylabel('Log Eff Moisture Permeability [log(kg/m^3)]')

  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()


def hysteresis(inv_alg_wp, inv_alg_dp,
               name='Hysteresis',
               Pc='inv_Pc',
               sat='inv_sat',
               seq='inv_seq',
               color=['k', 'b'],
               linestyle=['-','--'],
               marker=['o','^'],
               x1scale='log',
               y1scale=None,
               x2scale='log',
               y2scale='log',
               legend=['adsorption', 'desorption']):
  r"""
  Plot a montage of key saturation plots:
  - Sat(pc)
  - Sat(sim_step)
  - Pc(sat)
  - Pc(sim_step)

  Parameters
  ----------
  inv_alg : OpenPNM Algorithm Object
    The invasion algorithm for which the graphs are desired

  Examples
  --------

  >>> import OpenPNM
  >>> pn = OpenPNM.Network.TestNet()
  >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
  >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
  >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
  >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
  >>> IP = OpenPNM.Algorithms.InvasionPercolation(network=pn, name='IP')
  >>> IP.run(invading_phase=phase1, defending_phase=phase2, inlets=[pn.pores('top')], outlets=pn.pores('bottom'))
       IP algorithm at 0 % completion at 0.0 seconds
       IP algorithm at 20 % completion at 0.0 seconds
       IP algorithm at 40 % completion at 0.0 seconds
       IP algorithm at 60 % completion at 0.0 seconds
       IP algorithm at 100% completion at  0.0  seconds
  >>> OpenPNM.Postprocessing.Plots.drainage_curves(IP,timing='inv_time')

  """
  fig = _plt.figure(name, figsize=(6, 10), dpi=80,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

  pc_wp = inv_alg_wp._Pc
  w_wp = inv_alg_wp._moisturecontent
  try:
      eff_permeability_wp = inv_alg_wp._permeability
      _sp.log10(eff_permeability_wp)
  except:
      eff_permeability_wp = inv_alg_wp._conductance

  pc_dp = inv_alg_dp._Pc
  w_dp = inv_alg_dp._moisturecontent
  try:
      eff_permeability_dp = inv_alg_dp._permeability
      _sp.log10(eff_permeability_dp)
  except:
      eff_permeability_dp = inv_alg_dp._conductance

  p11, = ax1.plot(-pc_wp, w_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p12, = ax1.plot(-pc_dp, w_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax1.set_xlabel('Capillary Pressure [Pa]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax1.legend([p11, p12], legend)
  if x1scale is not None:
      ax1.set_xscale(x1scale)
      ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Log Moisture Content [log(kg/m^3)]')

  p21, = ax2.plot(-pc_wp, eff_permeability_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p22, = ax2.plot(-pc_dp, eff_permeability_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Effective Moisture Permeability [kg/msPa]')
  ax2.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax2.legend([p21, p22], legend)
  if x2scale is not None:
      ax2.set_xscale(x2scale)
      ax2.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y2scale is not None:
      ax2.set_yscale(y2scale)
      ax2.set_ylabel('Log Eff Moisture Permeability [log(kg/msPa)]')

  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()

def hysteresis_from_dict(data,
                         name='Hysteresis',
                         Pc='inv_Pc',
                         sat='inv_sat',
                         seq='inv_seq',
                         color=['k', 'b'],
                         linestyle=['-','--'],
                         marker=['o','^'],
                         x1scale='log',
                         y1scale=None,
                         x2scale='log',
                         y2scale='log',
                         legend=['adsorption', 'desorption']):
  r"""
  Plot a montage of key saturation plots:
  - Sat(pc)
  - Sat(sim_step)
  - Pc(sat)
  - Pc(sim_step)

  Parameters
  ----------
  data : a dictionary with keys = ['wetting', 'drying_wetting'].
    Each data[key] has the hygric property terms as its keys.

  Examples
  --------

  """
  ## EXTRACT data
  for key in data.keys():
      if 'drying' in key.split('_'):
          pc_dp = _sp.array(data[key]['pc'])
          w_dp = data[key]['moisture content']
          try:
              k_dp = data[key]['k_moisture']
          except:
              k_dp = data[key][0]
      else:
          pc_wp = _sp.array(data[key]['pc'])
          w_wp = data[key]['moisture content']
          try:
              k_wp = data[key]['k_moisture']
          except:
              k_wp = data[key][0]


  fig = _plt.figure(name, figsize=(6, 10), dpi=300,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

  p11, = ax1.plot(-pc_wp, w_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p12, = ax1.plot(-pc_dp, w_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax1.set_xlabel('Capillary Pressure [Pa]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax1.legend([p11, p12], legend)
  if x1scale is not None:
      ax1.set_xscale(x1scale)
      ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Log Moisture Content [log(kg/m^3)]')

  p21, = ax2.plot(-pc_wp, k_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p22, = ax2.plot(-pc_dp, k_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Effective Moisture Permeability [kg/msPa]')
  ax2.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax2.legend([p21, p22], legend)
  if x2scale is not None:
      ax2.set_xscale(x2scale)
      ax2.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y2scale is not None:
      ax2.set_yscale(y2scale)
      ax2.set_ylabel('Log Eff Moisture Permeability [log(kg/msPa)]')

  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()

def hysteresis_relative(inv_alg_wp, inv_alg_dp,
               name='Hysteresis',
               Pc='inv_Pc',
               sat='inv_sat',
               seq='inv_seq',
               color=['k', 'b'],
               linestyle=['-','--'],
               marker=['o','^'],
               x1scale='log',
               y1scale=None,
               x2scale='log',
               y2scale='log',
               legend=['adsorption', 'desorption']):
  r"""

  """
  fig = _plt.figure(name, figsize=(6, 10), dpi=80,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

  pc_wp = inv_alg_wp._Pc
  w_wp = inv_alg_wp._moisturecontent

  pc_dp = inv_alg_dp._Pc
  w_dp = inv_alg_dp._moisturecontent

  try:
      eff_permeability_wp = inv_alg_wp._corr_permeability
      eff_permeability_dp = inv_alg_dp._corr_permeability
  except:
      eff_permeability_wp = inv_alg_wp._rel_permeability
      eff_permeability_dp = inv_alg_dp._rel_permeability

  p11, = ax1.plot(-pc_wp, w_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p12, = ax1.plot(-pc_dp, w_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax1.set_xlabel('Capillary Pressure [Pa]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax1.legend([p11, p12], legend)
  if x1scale is not None:
      ax1.set_xscale(x1scale)
      ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Log Moisture Content [log(kg/m^3)]')

  p21, = ax2.plot(-pc_wp, eff_permeability_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  p22, = ax2.plot(-pc_dp, eff_permeability_dp, color=color[1],
                  linestyle=linestyle[1], marker=marker[1])
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Relative Moisture Permeability [-]')
  ax2.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax2.legend([p21, p22], legend)
  if x2scale is not None:
      ax2.set_xscale(x2scale)
      ax2.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y2scale is not None:
      ax2.set_yscale(y2scale)
      ax2.set_ylabel('Log Rel Moisture Permeability [-]')

  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()

def _extract_pt(inv_alg, Pc='inv_Pc', sat='inv_sat', seq='inv_seq'):
    r"""
    Extract whole (for pores and throats -> pt) pc, sat, and sequences
    """
    inv_throats = inv_alg.toindices(inv_alg['throat.'+seq]>0)
    sort_seq = _sp.argsort(inv_alg['throat.'+seq][inv_throats])
    inv_throats = inv_throats[sort_seq]
    tpc_unik = _sp.unique(inv_alg['throat.'+Pc][inv_throats])
    tsat_unik = _sp.unique(inv_alg['throat.'+sat][inv_throats])
    tseq_unik = _sp.unique(inv_alg['throat.'+seq][inv_throats])
    t_unik = _sp.array(list(zip(tpc_unik, tsat_unik, tseq_unik)))

    inv_pores = inv_alg.toindices(inv_alg['pore.'+seq]>0)
    sort_seq = _sp.argsort(inv_alg['pore.'+seq][inv_pores])
    inv_pores = inv_pores[sort_seq]
    ppc_unik = _sp.unique(inv_alg['pore.'+Pc][inv_pores])
    psat_unik = _sp.unique(inv_alg['pore.'+sat][inv_pores])
    pseq_unik = _sp.unique(inv_alg['pore.'+seq][inv_pores])
    p_unik = _sp.array(list(zip(ppc_unik, psat_unik, pseq_unik)))
    unik = _sp.concatenate((t_unik, p_unik))

    unik = unik[unik[:,1].argsort()]

    pt_sat = _sp.unique(unik[:,1])
    pt_seq = _sp.unique(unik[:,2])
    pt_pc = _sp.zeros_like(pt_sat)
    j = 0
    nn = _sp.shape(unik)[0]
    for i in range(0,nn-1):
        if unik[i,2] == unik[i+1,2]:
            pt_pc[j] = max(unik[i,0], unik[i+1,0])
            j += 1
        elif _sp.mod(i, 2)==0:
            pt_pc[j] = unik[i,0]
            j += 1
        elif unik[i,2] != unik[i-1,2]:
            pt_pc[j] = unik[i,0]
            j += 1
    if i+1 == nn-1:
        pt_pc[j] = unik[i+1,0]

    return (pt_seq, pt_pc, pt_sat)

def plot_wp(inv_alg_wp,
               name='Wetting',
               Pc='inv_Pc',
               sat='inv_sat',
               seq='inv_seq',
               color=['k'],
               linestyle=['-'],
               marker=['o'],
               x1scale='log',
               y1scale=None,
               x2scale='log',
               y2scale='log',
               legend=None):
  r"""
  Plot the moisture retention curve and permeability of wetting percolation algorithm in one figure of two subplots
  """
  if legend is None:
      if name.lower() == 'wetting':
          legend = ['adsorption']
      elif name.lower() == 'imbibition':
          legend = ['imbibition']

  fig = _plt.figure(name, figsize=(6, 10), dpi=80,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(211)   # top
  ax2 = fig.add_subplot(212)   # bottom

  pc_wp = inv_alg_wp._Pc
  w_wp = inv_alg_wp._moisturecontent
  try:
      eff_permeability_wp = inv_alg_wp._permeability
      _sp.log10(eff_permeability_wp)
  except:
      eff_permeability_wp = inv_alg_wp.result['k_moisture']
#      eff_permeability_wp = inv_alg_wp._conductance


  p11, = ax1.plot(-pc_wp, w_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  ax1.set_xlabel('Capillary Pressure [Pa]')
  ax1.set_ylabel('Moisture Content [kg/m^3]')
  ax1.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax1.legend([p11], legend)
  if x1scale is not None:
      ax1.set_xscale(x1scale)
      ax1.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Log Moisture Content [log(kg/m^3)]')

  p21, = ax2.plot(-pc_wp, eff_permeability_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  ax2.set_xlabel('Capillary Pressure [Pa]')
  ax2.set_ylabel('Effective Moisture Permeability [kg/msPa]')
  ax2.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax2.legend([p21], legend, loc=1)
  if x2scale is not None:
      ax2.set_xscale(x2scale)
      ax2.set_xlabel('Log Capillary Pressure [log(-Pa)]')
  if y2scale is not None:
      ax2.set_yscale(y2scale)
      ax2.set_ylabel('Log Eff Moisture Permeability [log(kg/msPa)]')

  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)
  ax2.grid(True)

  fig.show()

def k_w(inv_alg_wp,
        name='Permeability(moisture)',
        color=['k'],
        linestyle=['-'],
        marker=['o'],
        y1scale=None,
        legend=['adsorption']):
  r"""
  Plot the permeability w.r.t the moisture content
  """
  fig = _plt.figure(name, figsize=(6, 10), dpi=80,
                    facecolor='w', edgecolor='k')
  ax1 = fig.add_subplot(111)   # top

  pc_wp = inv_alg_wp._Pc
  w_wp = inv_alg_wp._moisturecontent
  try:
      eff_permeability_wp = inv_alg_wp._permeability
      _sp.log10(eff_permeability_wp)
  except:
      eff_permeability_wp = inv_alg_wp._conductance


  p11, = ax1.plot(w_wp, eff_permeability_wp, color=color[0],
                  linestyle=linestyle[0], marker=marker[0])
  ax1.set_xlabel('Moisture content [kg/m^3]')
  ax1.set_ylabel('Permeability [log(kg/msPa)]')
  ax1.set_xlim([1.01*min(-pc_wp), 0.99*max(-pc_wp)])
  ax1.legend([p11], legend)
  if y1scale is not None:
      ax1.set_yscale(y1scale)
      ax1.set_ylabel('Permeability [log(kg/msPa)]')
  fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
  ax1.grid(True)

  fig.show()

def plot_2scales(x, y1, y2,
                 title='w and k',
                 xlabel='-Capillary Pressure [-Pa]',
                 y1label='Moisture Content [kg/m^3]',
                 y2label='Log Capillary Pressure [log(-Pa)]',
                 xscale='log', y1scale=None, y2scale='log',
                 color=['k', 'b'], linestyle=['-','-'], marker=['o','o']):
  r"""
  Plot a figure with 2 graphs and 2 y-scales

  Arguments:
  ----------
  x:    -pc (array-like)
  y1:   moisture content (array-like)
  y2:   permeability (array-like)
  """
  fig, ax1 = _plt.subplots()

  ax1.plot(x, y1, color=color[0], linestyle=linestyle[0], marker=marker[0])
  ax1.set_xlabel(xlabel)
  # Make the y-axis label and tick labels match the line color.
  ax1.set_ylabel(y1label, color=color[0])
  ax1.set_xscale(xscale)
  if y1scale is not None:
      ax1.set_yscale(y1scale)
  for tl in ax1.get_yticklabels():
      tl.set_color(color[0])

  ax2 = ax1.twinx()

  ax2.plot(x, y2, color=color[1], linestyle=linestyle[1], marker=marker[1])
  ax2.set_ylabel(y2label, color=color[1])
  ax2.set_xscale(xscale)
  ax2.set_yscale(y2scale)
  for tl in ax2.get_yticklabels():
      tl.set_color(color[1])

  ax1.grid(True)
  ax2.grid(True)
  _plt.show()


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)

