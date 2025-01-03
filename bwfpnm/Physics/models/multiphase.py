r"""
===============================================================================
Submodule -- diffusive_conductance
===============================================================================

"""
import numpy as np
import scipy as sp
import re


def mixed_conductance(physics,
                      phase,
                      network,
                      diffusion=True,
                      pore_occupancy='pore.occupancy',
                      throat_occupancy='throat.occupancy',
                      thydraulic_conductance='throat.hydraulic_conductance',
                      tdiffusive_conductance='throat.diffusive_conductance',
                      **kwargs):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.

    This method "closes" conduits that are not sufficiently filled with the
    specified phase by multiplying the original conductance by a very small *factor*.

    phase: MOISTURE object

    """
    phydraulic_conductance = thydraulic_conductance + '_pore'
    pdiffusive_conductance = tdiffusive_conductance + '_pore'

    _physics = physics._net._physics
    for phys in _physics:
        words = phys.name.lower().split('_')
        words += re.findall('[A-Z][^A-Z]*', phys.name)
        words = [word.lower() for word in words]
        if 'water' in words:
            phys_wtr = phys
        elif 'vapour' in words:
            phys_vpr = phys
        elif 'vap' in words:
            phys_vpr = phys

    _phases = physics._net._phases
    for _phase in _phases:
        if 'water' == _phase.name:
            liquid = _phase
            break

    p1 = network['throat.conns'][:,0]
    p2 = network['throat.conns'][:,1]
    # wet pore = True, dry pore = False
    t_liquid = sp.array(liquid[throat_occupancy], dtype=bool)
    p1_liquid = sp.array(liquid[pore_occupancy][p1], dtype=bool)
    p2_liquid = sp.array(liquid[pore_occupancy][p2], dtype=bool)

    # get all conductances from each physics of phases
    gp = phys_wtr[phydraulic_conductance]     # pores
    gt = phys_wtr[thydraulic_conductance]     # throats
    gp_v = phys_vpr[pdiffusive_conductance]
    gt_v = phys_vpr[tdiffusive_conductance]

    # select appropriate (liquid/vapour) conductances based on their occupancy
    gp1 = gp[:,0]*p1_liquid + gp_v[:,0]*(~p1_liquid)
    gp2 = gp[:,1]*p2_liquid + gp_v[:,1]*(~p2_liquid)
    gt = gt*t_liquid + gt_v*(~t_liquid)
#    value = (1 / gt + 1 / gp1 + 1 / gp2) ** (-1)
#    print(gp1, gp2, gt, value)
    # For throat/pore network: volumeless element gives no restriction to flow
    # set conductance = infinity
#    gp1[gp1 <= 1e-100] = sp.inf
#    gp2[gp2 <= 1e-100] = sp.inf
#    gt[gt <= 1e-100] = sp.inf

    # Total g: harmonic mean: value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = harmonic_mean([gt,gp1,gp2])

    if not diffusion:
        conduits = 1*t_liquid + 1*p1_liquid + 1*p2_liquid
        value[conduits<3] = 0

    return value

def mixed_surf_conduc(physics,
                      phase,
                      network,
                      diffusion=True,
                      pore_occupancy='pore.occupancy',
                      throat_occupancy='throat.occupancy',
                      thydraulic_conductance='throat.hydraulic_conductance',
                      tdiffusive_conductance='throat.diffusive_conductance',
                      tsurf_diff_cond='throat.surface_conductance',
                      **kwargs):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.

    This method "closes" conduits that are not sufficiently filled with the
    specified phase by multiplying the original conductance by a very small *factor*.

    phase: MOISTURE object

    """
    diffusion = True
    phydraulic_conductance = thydraulic_conductance + '_pore'
    pdiffusive_conductance = tdiffusive_conductance + '_pore'
    psurf_diff_cond = tsurf_diff_cond + '_pore'

    _physics = physics._net._physics
    for phys in _physics:
        words = phys.name.lower().split('_')
        words += re.findall('[A-Z][^A-Z]*', phys.name)
        words = [word.lower() for word in words]
        if 'water' in words:
            phys_wtr = phys
        elif 'vapour' in words:
            phys_vpr = phys
        elif 'vap' in words:
            phys_vpr = phys

    _phases = physics._net._phases
    for _phase in _phases:
        if 'water' in _phase.name.split('_'):
            liquid = _phase
            break

    p1 = network['throat.conns'][:,0]
    p2 = network['throat.conns'][:,1]
    # wet pore = True, dry pore = False
    t_liquid = sp.array(liquid[throat_occupancy], dtype=bool)
    p1_liquid = sp.array(liquid[pore_occupancy][p1], dtype=bool)
    p2_liquid = sp.array(liquid[pore_occupancy][p2], dtype=bool)
    conduits = 1*t_liquid + 1*p1_liquid + 1*p2_liquid

    # get all conductances from each physics of phases
    gl_p = phys_wtr[phydraulic_conductance]     # pores
    gl_t = phys_wtr[thydraulic_conductance]     # throats
    gv_p = phys_vpr[pdiffusive_conductance]
    gv_t = phys_vpr[tdiffusive_conductance]
    glsurf_p = phys_wtr[psurf_diff_cond]
    glsurf_t = phys_wtr[tsurf_diff_cond]

    # copy the hydraulic conductances of wet pores.gt[drypores] = 0
    gp1 = gl_p[:,0]*p1_liquid               # hydraulic conductance of wet pore
    gp1_v = gv_p[:,0]*(~p1_liquid)          # diffusive conductance of dry pore
    gp1_surf = glsurf_p[:,0]*(~p1_liquid)   # surface conductance of dry pore
    gp2 = gl_p[:,1]*p2_liquid
    gp2_v = gv_p[:,1]*(~p2_liquid)
    gp2_surf = glsurf_p[:,1]*(~p2_liquid)
    gt = gl_t*t_liquid
    gt_v = gv_t*(~t_liquid)
    gt_surf = glsurf_t*(~t_liquid)

    # treating dry pores (with gv and glsurf)
    val_v, val_surf = sp.zeros_like(gt), sp.zeros_like(gt)
    val = sp.zeros_like(gt)
    # condition 0: all dry elements
    cond0 = conduits == 0 #
    if sp.any(cond0):
        val_v[cond0] = harmonic_mean([gp1_v[cond0], gp2_v[cond0], gt_v[cond0]])
        val_surf[cond0] = harmonic_mean([gp1_surf[cond0], gp2_surf[cond0],
                                         gt_surf[cond0]])
        val[cond0] = val_v[cond0]+val_surf[cond0]
    # condition 1: 1 wet element
    cond11 = (conduits == 1) * (gp1 != 0)
    if sp.any(cond11):
        val_v[cond11] = harmonic_mean([gp2_v[cond11], gt_v[cond11]])
        val_surf[cond11] = harmonic_mean([gp2_surf[cond11], gt_surf[cond11]])
        val[cond11] = val_v[cond11] + val_surf[cond11]
        val[cond11] = harmonic_mean([val[cond11], gp1[cond11]])

    cond12 = (conduits == 1) * (gp2 != 0)
    if sp.any(cond12):
        val_v[cond12] = harmonic_mean([gp1_v[cond12], gt_v[cond12]])
        val_surf[cond12] = harmonic_mean([gp1_surf[cond12], gt_surf[cond12]])
        val[cond12] = val_v[cond12] + val_surf[cond12]
        val[cond12] = harmonic_mean([val[cond12], gp2[cond12]])

    cond13 = (conduits == 1) * (gt != 0)
    if sp.any(cond13):
        val[cond13] = harmonic_mean([gt[cond13], gp1_v[cond13] + gp1_surf[cond13],
                                   gp2_v[cond13] + gp2_surf[cond13]])
    # condition 2: 2 wet elements
    cond21 = (conduits == 2) * (gp1 == 0)
    if sp.any(cond21):
        val[cond21] = harmonic_mean([gt[cond21], gp2[cond21],
                                   gp1_v[cond21] + gp1_surf[cond21]])
    cond22 = (conduits == 2) * (gp2 == 0)
    if sp.any(cond22):
        val[cond22] = harmonic_mean([gt[cond22], gp1[cond22],
                                   gp2_v[cond22] + gp2_surf[cond22]])
    cond23 = (conduits == 2) * (gt == 0)
    if sp.any(cond23):
        val[cond23] = harmonic_mean([gp1[cond23], gp2[cond23],
                                   gt_v[cond23] + gt_surf[cond23]])
    # condition 3: all wet elements
    cond3 = conduits == 3
    if sp.any(cond3):
        val[cond3] = harmonic_mean([gp1[cond3], gp2[cond3], gt[cond3]])

    cond = 1*cond0 + 1*cond11 + 1*cond12 + 1*cond13 + 1*cond21 + 1*cond22 + 1*cond23 + 1*cond3
    valid = sp.all(cond==1)
    if not valid:
        print('Invalid operation in multiphase.mixed_surf_conduc' )
    # For throat/pore network: volumeless element gives no restriction to flow
    # set conductance = infinity
    if not diffusion:
        val[conduits<3] = 0

    return val

def harmonic_mean(array):
    array = sp.array(array)
    # Total g: harmonic mean: value = (1/gt + 1/gp1 + 1/gp2)**(-1)
    value = sp.sum(1.0/array, axis=0) # any zero in gt or gp will result in zero mean
    value = value**(-1)
    if sp.any(value==sp.inf):
        print('Error: infinity-valued moisture conduit conductance. Something goes wrong!')
    return value

def single_conductance_pore(physics,
                           network,
                           factor=1e-40,
                           throat_occupancy='throat.occupancy',
                           pore_occupancy='pore.occupancy',
                           pconductance='throat.hydraulic_conductance_pore',
                           tconductance='throat.hydraulic_conductance',
                           psurfcond='throat.surface_conductance_pore',
                           tsurfcond='throat.surface_conductance',
                           **kwargs):
    r"""
    Add a new multiphase conductance property to the conduits of network, where a
    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.

    This method "closes" conduits that are not sufficiently filled with the
    specified phase by multiplying the original conductance by a very small *factor*.

    """
    p1 = network['throat.conns'][:,0]
    p2 = network['throat.conns'][:,1]
    t_open = sp.array(physics[throat_occupancy], dtype=bool)
    p1_open = sp.array(physics[pore_occupancy][p1], dtype=bool)
    p2_open = sp.array(physics[pore_occupancy][p2], dtype=bool)

    gp = physics[pconductance]
    gt = physics[tconductance]*t_open
    gp1 = gp[:,0]*p1_open
    gp2 = gp[:,1]*p2_open

    try:    # surface coductances, only valid for phys_water
        gpsurf = physics[psurfcond]
        gtsurf = physics[tsurfcond]
        gt[~t_open] = gtsurf[~t_open]
        gp1[~p1_open] = gpsurf[:,0][~p1_open]
        gp2[~p2_open] = gpsurf[:,1][~p2_open]
        surflow = True
    except:
        surflow = False

#    # For throat/pore network: volumeless element gives no restriction to flow
#    # set conductance = infinity
#    gp1[gp1 <= 1e-100] = sp.inf
#    gp2[gp2 <= 1e-100] = sp.inf
#    gt[gt <= 1e-100] = sp.inf

    # Total g: harmonic mean (sum of resistances)
    value = harmonic_mean([gp1, gt, gp2])

    if not surflow and ('water' in physics.name.split('_')):
        open_conduits = p1_open & t_open & p2_open
        closed_conduits = ~open_conduits
#        span, dead, isolated = network.cluster_types(mask=open_conduits,
#                                                     save=False)
#        if len(span) > 0 and len(isolated) > 0:
#            # trimming the isolated conduits to avoid singular matrix, but
#            # the X is still nan even the spanning cluster exists.
#            trim_pores = sp.hstack(isolated)
#            t_ind = network.find_neighbor_throats(pores=trim_pores)
#            value[t_ind] = 0
#            # make value[throat[isolated]] = 0 or try use factor = 1e-40
        value[closed_conduits] = 0

#       value = value * open_conduits + value * closed_conduits * factor
    return value


# def mixed_conductance(physics,
#                      phase,
#                      network,
#                      throat_conductance='mixed',
#                      throat_occupancy='throat.occupancy',
#                      pore_occupancy='pore.occupancy',
#                      phydraulic_conductance='pore.hydraulic_conductance',
#                      pdiffusive_conductance='pore.diffusive_conductance',
#                      thydraulic_conductance='throat.hydraulic_conductance',
#                      tdiffusive_conductance='throat.diffusive_conductance',
#                      **kwargs):
#    r"""
#    Add a new multiphase conductance property to the conduits of network, where a
#    conduit is ( 1/2 pore - full throat - 1/2 pore ) based on the areas.
#
#    This method "closes" conduits that are not sufficiently filled with the
#    specified phase by multiplying the original conductance by a very small *factor*.
#
#    Parameters
#    ----------
#    network : OpenPNM Network Object
#
#    phase : OpenPNM Phase Object
#        The phase of interest
#
#    occupied_condition : 'occupancy'
#        The name of the pore and throat property that dictates whether conduit is "closed" or not
#
#    Notes
#    -----
#    This function requires that all the necessary phase properties already be
#    calculated.
#
#    """
#    _physics = physics._net._physics
#    for phys in _physics:
#        words = phys.name.lower().split('_')
#        words += re.findall('[A-Z][^A-Z]*', phys.name)
#        words = [word.lower() for word in words]
#        if 'water' in words:
#            liquid = phys
#        elif 'vapour' in words:
#            vapour = phys
#        elif 'vap' in words:
#            vapour = phys
##    throats = phase.throats(physics.name)
#    throats = network.throats()
#    connected_pores = network.find_connected_pores(throats)
#    pores_1 = connected_pores[:, 0]
#    pores_2 = connected_pores[:, 1]
#
#    if throat_conductance == 'mixed':
#        # get all conductances from each physics of phases
#        gp1 = liquid[phydraulic_conductance][pores_1]
#        gp1_v = vapour[pdiffusive_conductance][pores_1]
#        gp2 = liquid[phydraulic_conductance][pores_2]
#        gp2_v = vapour[pdiffusive_conductance][pores_2]
#        gt = liquid[thydraulic_conductance]
#        gt_v = vapour[tdiffusive_conductance]
#    else:
#        try:
#            gp1 = phase[throat_conductance][pores_1]
#            gp1_v = gp1
#            gp2 = phase[throat_conductance][pores_2]
#            gp2_v = gp2
#            gt = phase[throat_conductance]
#            gt_v = gt
#        except:
#            pass
#
#    # select appropriate conductances based on their occupancy
#    _phases = physics._net._phases
#    for _phase in _phases:
#        if 'water' == _phase.name:
#            liquid = _phase
#            break
#
#    throats_liquid = sp.array(liquid[throat_occupancy], dtype=bool)
#    pores_1_liquid = sp.array(liquid[pore_occupancy][pores_1], dtype=bool)
#    pores_2_liquid = sp.array(liquid[pore_occupancy][pores_2], dtype=bool)
#
#    gp1 = gp1*pores_1_liquid + gp1_v*(-pores_1_liquid)
#    gp2 = gp2*pores_2_liquid + gp2_v*(-pores_2_liquid)
#    gt = gt*throats_liquid + gt_v*(-throats_liquid)
#
#    # For volumeless pores/throats set conductance = infinity
#    gp1[gp1 <= 1e-100] = sp.inf
#    gp2[gp2 <= 1e-100] = sp.inf
#    gt[gt <= 1e-100] = sp.inf
#
#    # Total g
#    value = (1/gt + 1/gp1 + 1/gp2)**(-1)
##    print('pconductance: ', sp.amin(sp.r_[gp1, gp2]))
##    print('tconductance: ', gt[0:10])
#    return value
