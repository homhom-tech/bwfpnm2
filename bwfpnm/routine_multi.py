# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 10:09:46 2015

@author: islah
"""

import bwfpnm as bpnm
import scipy as sp
import copy

ctrl = bpnm.Base.Controller()
ctrl.loglevel = 40


def net_geo(netnames=None, geonames=None, shapes=None, z=None):
    r'''
    all arguments must be dictionary with two keys: 'macro' and 'micro'
    '''
    NetClass = bpnm.Network.RegularLattice
    GeoClass = bpnm.Geometry.Cylinder_Sphere
    pn = {'macro': 0, 'micro': 0}
    geo = copy.deepcopy(pn)
    if netnames is None:
        netnames = {'macro': 'NetMacro', 'micro': 'NetMicro'}
    if geonames is None:
        geonames = {'macro': 'GeoMacro', 'micro': 'GeoMicro'}
    if shapes is None:
        shapes = {'macro': [2, 2, 2], 'micro': [2, 2, 2]}
    if z is None:
        z = {'macro': 6, 'micro': 6}

    for key in pn.keys():
        pn[key] = NetClass(name=netnames[key], shape=shapes[key],
                           connectivity=z[key])
        if key == 'micro':
            pn[key].add_boundaries(labels=['left', 'right'])
        Ps = pn[key].pores()
        Ts = pn[key].throats()
        geo[key] = GeoClass(network=pn[key], pores=Ps, throats=Ts,
                            name=geonames[key], scale=key)
    return (pn, geo)


def phase(pn, phasenames=None):
    r'''
    all arguments must be dictionary with two keys: 'macro' and 'micro'
    '''
    WaterClass = bpnm.Phases.Water
    VapourClass = bpnm.Phases.Vapour
    MoistureClass = bpnm.Phases.GenericPhase
    if phasenames is None:
        phasenames = {'macro': 'Macro', 'micro': 'Micro'}
    water = {'macro': 0, 'micro': 0}
    vapour = copy.deepcopy(water)
    moisture = copy.deepcopy(water)

    for key in water.keys():
        water[key] = WaterClass(name='Water'+phasenames[key], network=pn[key])
        vapour[key] = VapourClass(name='Vap'+phasenames[key], network=pn[key])
        moisture[key] = MoistureClass(name='Moist'+phasenames[key],
                                      network=pn[key])
        vapour[key]['pore.water_density'] = water[key]['pore.density']
        moisture[key]['pore.temperature'] = 293.15

    return (water, vapour, moisture)


def phys_phase(pn, water, vapour, moisture, param='radius'):
    r'''
    all arguments must be dictionary with two keys: 'macro' and 'micro'
    '''
    PhysWaterClass = bpnm.Physics.Standard_Topology
    PhysVapClass = bpnm.Physics.GenericPhysics
    PhysMoistClass = bpnm.Physics.GenericPhysics

    names = {'macro': 'Macro', 'micro': 'Micro'}
    phys_water = {'macro': 0, 'micro': 0}
    phys_vapour = {'macro': 0, 'micro': 0}
    phys_moisture = {'macro': 0, 'micro': 0}
    for key in phys_water.keys():
        ps = pn[key].pores()
        ts = pn[key].throats()
        phys_water[key] = PhysWaterClass(network=pn[key], phase=water[key],
                                         pores=ps, throats=ts, param=param,
                                         dynamic_data=True,
                                         name='Phys_Water_'+names[key])
        phys_vapour[key] = PhysVapClass(network=pn[key], phase=vapour[key],
                                        pores=ps, throats=ts,
                                        dynamic_data=True,
                                        name='Phys_Vap_'+names[key])
        phys_moisture[key] = PhysMoistClass(network=pn[key],
                                            phase=moisture[key],
                                            pores=ps, throats=ts,
                                            dynamic_data=True,
                                            name='Phys_Moist_'+names[key])
    return (phys_water, phys_vapour, phys_moisture)


def algorithm(pn, water, vapour, inv_points=None, npts=None, drying=True):
    r'''
    pn, water, and vapour arguments must be dictionary with
    two keys: 'macro' and 'micro'

    return (alg_wp, alg_dp) if drying is True, else return alg_wp
    '''
    AdsorpClass = bpnm.Algorithms.WettingPercolation
    DesorpClass = bpnm.Algorithms.DryingPercolation
    if npts is None:
        npts = 10

    names = {'macro': 'Macro', 'micro': 'Micro'}
    alg_wp = {'macro': 0, 'micro': 0}
    alg_dp = {'macro': 0, 'micro': 0}
    for key in alg_wp.keys():
        alg_wp[key] = AdsorpClass(network=pn[key],
                                  invading_phase=water[key],
                                  defending_phase=vapour[key],
                                  name='WettingPercolation'+names[key])
        inv_sites = pn[key]['pore.inlet'] + pn[key]['pore.outlet']    # bool
        if inv_points is None:
            pc_list = sp.r_[water[key]['pore.capillary_pressure'],
                            water[key]['throat.capillary_pressure']]
            pc_list.sort()
        alg_wp[key].run(inlets=None, npts=npts, inv_points=pc_list)

        if drying:
            alg_dp[key] = DesorpClass(network=pn[key],
                                      invading_phase=vapour[key],
                                      defending_phase=water[key],
                                      name='DryingPercolation'+names[key])
            alg_dp[key].run(inlets=inv_sites, npts=npts, inv_points=pc_list)

    if drying:
        return (alg_wp, alg_dp)
    else:
        return alg_wp


def permeability(pn, alg_wp, alg_dp, water, vapour, moisture,
                 phys_vapour, phys_moisture, w_sat,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc =1,
                 **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    names = {'macro': 'Macro', 'micro': 'Micro'}
    alg_flow_wp = {'macro': 0, 'micro': 0, 'multi':0}
    alg_flow_dp = {'macro': 0, 'micro': 0, 'multi':0}
    keys = ['micro', 'macro'] #, 'multi']       # to ensure started from micro
    key0 = None
    for key in keys:
        print(key+'-scale network is calculated')
        sat_wp, sat_dp = [], []
        sat_wp_surf, sat_dp_surf = [], []
        sat_wp_moist, sat_dp_moist = [], []
        max_norm_res = 0
        eff_perm_moisture_wp = {'0': [], '1': [], '2': []}
        eff_perm_moisture_dp = {'0': [], '1': [], '2': []}

        if key == 'multi':
            key = 'macro'
            key0 = 'multi'

        p_volumes = pn[key]['pore.volume']
        t_volumes = pn[key]['throat.volume']
        volume_total = sum(p_volumes) + sum(t_volumes)
        if key0 == 'multi':
            Nt = pn[key].num_throats()
            p_volumes += Nt*pn['micro']['pore.volume']
            t_volumes += Nt*pn['micro']['throat.volume']
            volume_total = sum(p_volumes) + sum(t_volumes)

        lr = sp.arange(-10, -1)
        r = sp.power(10, lr)
        pc = -2*water[key]['pore.surface_tension'][0]/r

        Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
        Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
        Pc_wp = sp.around(Pc_wp, 3)
        Pc_dp = Pc_wp[::-1]

    #    Pc_min, Pc_max = sp.amin(alg_wp._inv_points), sp.amax(alg_wp._inv_points)
    #    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
    #    Pc_wp = sp.around(Pc_wp, 3)
    #    Pc_dp = Pc_wp[::-1]

        for Pc_step_wp, Pc_step_dp in list(zip(Pc_wp, Pc_dp)):
            alg_wp[key].return_results(Pc=Pc_step_wp, occupancy='occupancy_wp')
            alg_dp[key].return_results(Pc=Pc_step_dp, occupancy='occupancy_dp')

            p_occ_wp = water[key]['pore.occupancy_wp']
            t_occ_wp = water[key]['throat.occupancy_wp']
            p_occ_dp = water[key]['pore.occupancy_dp']
            t_occ_dp = water[key]['throat.occupancy_dp']

            if key0 == 'multi':
                p_occ_wp += Nt*water['micro']['pore.occupancy_wp']
                t_occ_wp += Nt*water['micro']['throat.occupancy_wp']
                p_occ_dp += Nt*water['micro']['pore.occupancy_dp']
                t_occ_dp += Nt*water['micro']['throat.occupancy_dp']

            volume_p_wp = sum(p_occ_wp*p_volumes)
            volume_t_wp = sum(t_occ_wp*t_volumes)
            volume_p_dp = sum(p_occ_dp*p_volumes)
            volume_t_dp = sum(t_occ_dp*t_volumes)

            saturation_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            saturation_dp = (volume_p_dp + volume_t_dp)/(volume_total)
            sat_wp.append(saturation_wp)
            sat_dp.append(saturation_dp)
            if printstatus:
                print('WP_saturation: %.3f' % saturation_wp,
                      '\t Pc: %.3f' %Pc_step_wp)
                print('DP_saturation: %.3f' % saturation_dp,
                      '\t Pc: %.3f' %Pc_step_dp)
                print('WP_volume: ', sum(p_occ_wp),
                      '\t throat: ', sum(t_occ_wp))
                print('DP_volume: ', sum(p_occ_dp),
                      '\t throat: ', sum(t_occ_dp))
                print('WP_water occupancy: ', p_occ_wp,
                      '\t throat: ', t_occ_wp)
                print('DP_water occupancy: ', p_occ_dp,
                      '\t throat: ', t_occ_dp, '\n')

            # Surface adsorption: p/t.surface_thickness & surface_volume
            if surface_ad:
                phys_vapour[key].models.add(propname='throat.surface_thickness_wp',
                                       model=pm.surface_adsorption.tstat_thickness,
                                       pc=Pc_step_wp)
                phys_vapour[key].models.add(propname='pore.surface_thickness_wp',
                                       model=pm.surface_adsorption.pstat_thickness,
                                       pc=Pc_step_wp)
                phys_vapour[key].models.add(propname='throat.surface_thickness_dp',
                                       model=pm.surface_adsorption.tstat_thickness,
                                       pc=Pc_step_dp,
                                       throat_occupancy='throat.occupancy_dp')
                phys_vapour[key].models.add(propname='pore.surface_thickness_dp',
                                       model=pm.surface_adsorption.pstat_thickness,
                                       pc=Pc_step_dp,
                                       pore_occupancy='pore.occupancy_dp')

                phys_vapour[key].models.add(propname='throat.surface_volume_wp',
                                       model=pm.surface_adsorption.tvolume,
                                       film_thickness='throat.surface_thickness_wp')
                phys_vapour[key].models.add(propname='pore.surface_volume_wp',
                                       model=pm.surface_adsorption.pvolume,
                                       film_thickness='pore.surface_thickness_wp')
                phys_vapour[key].models.add(propname='throat.surface_volume_dp',
                                       model=pm.surface_adsorption.tvolume,
                                       throat_occupancy='throat.occupancy_dp',
                                       film_thickness='throat.surface_thickness_dp')
                phys_vapour[key].models.add(propname='pore.surface_volume_dp',
                                       model=pm.surface_adsorption.pvolume,
                                       pore_occupancy='pore.occupancy_dp',
                                       film_thickness='pore.surface_thickness_dp')

                volume_p_wp = sum(phys_vapour[key]['pore.surface_volume_wp'])
                volume_t_wp = sum(phys_vapour[key]['throat.surface_volume_wp'])
                volume_p_dp = sum(phys_vapour[key]['pore.surface_volume_dp'])
                volume_t_dp = sum(phys_vapour[key]['throat.surface_volume_dp'])

                sat_surf_wp = (volume_p_wp + volume_t_wp)/(volume_total)
                sat_surf_dp = (volume_p_dp + volume_t_dp)/(volume_total)
                sat_wp_surf.append(sat_surf_wp)
                sat_dp_surf.append(sat_surf_dp)

                sat_wp[-1] += sat_surf_wp   # update total saturations
                sat_dp[-1] += sat_surf_dp

                print('pthickness wp: ',
                      phys_vapour[key]['pore.surface_thickness_wp'][0])
                print('pradius wp: ',
                      pn[key]['pore.diameter'][0]/2)
                print('pvol surf wp: ',
                      phys_vapour[key]['pore.surface_volume_wp'][0])
                print('pvol wp: ',
                      p_volumes[0])
                print('psat wp: ',
                      phys_vapour[key]['pore.surface_volume_wp'][0]/p_volumes[0])
    #            print('tthickness wp: ',
    #                  phys_vapour[key]['throat.surface_thickness_wp'])
    #            print('tsat wp: ',
    #                  phys_vapour[key]['throat.surface_volume_wp']/t_volumes)

            if moist_volume:
                phys_vapour[key].models.add(propname='throat.moist_volume_wp',
                                       model=pm.volume_moisture.tvolume,
                                       pc=Pc_step_wp)
                phys_vapour[key].models.add(propname='pore.moist_volume_wp',
                                       model=pm.volume_moisture.pvolume,
                                       pc=Pc_step_wp)
                phys_vapour[key].models.add(propname='throat.moist_volume_dp',
                                       model=pm.volume_moisture.tvolume,
                                       pc=Pc_step_dp,
                                       throat_occupancy='throat.occupancy_dp')
                phys_vapour[key].models.add(propname='pore.moist_volume_dp',
                                       model=pm.volume_moisture.pvolume,
                                       pc=Pc_step_dp,
                                       pore_occupancy='pore.occupancy_dp')

                volume_p_wp = sum(phys_vapour[key]['pore.moist_volume_wp'])
                volume_t_wp = sum(phys_vapour[key]['throat.moist_volume_wp'])
                volume_p_dp = sum(phys_vapour[key]['pore.moist_volume_dp'])
                volume_t_dp = sum(phys_vapour[key]['throat.moist_volume_dp'])

                sat_moist_wp = (volume_p_wp + volume_t_wp)/(volume_total)
                sat_moist_dp = (volume_p_dp + volume_t_dp)/(volume_total)
                sat_wp_moist.append(sat_moist_wp)
                sat_dp_moist.append(sat_moist_dp)

                sat_wp[-1] += sat_moist_wp   # update total saturations
                sat_dp[-1] += sat_moist_dp

                print('moist vol: ', phys_vapour[key]['throat.moist_volume_wp'])
                print('moist sat: ',
                      phys_vapour[key]['throat.moist_volume_wp']/t_volumes)

            # Update vapour permeability for all pores & throats for wp & dp
            phys_vapour[key].models.add(propname='throat.diffusive_conductance_wp',
                                   model=pm.diffusive_conductance.tbulk_diffusion,
                                   pc=Pc_step_wp, knudsen=knudsen)
            phys_vapour[key].models.add(propname='pore.diffusive_conductance_wp',
                                   model=pm.diffusive_conductance.pbulk_diffusion,
                                   pc=Pc_step_wp, knudsen=knudsen)
            phys_vapour[key].models.add(propname='throat.diffusive_conductance_dp',
                                   model=pm.diffusive_conductance.tbulk_diffusion,
                                   pc=Pc_step_dp, knudsen=knudsen)
            phys_vapour[key].models.add(propname='pore.diffusive_conductance_dp',
                                   model=pm.diffusive_conductance.pbulk_diffusion,
                                   pc=Pc_step_dp, knudsen=knudsen)
            phys_vapour[key].models.regenerate()
            # Calculate conduit conductances as a function of water distribution
            phys_moisture[key].models.add(propname='throat.conduit_conductance_wp',
                                     model=pm.multiphase.mixed_conductance,
                                     throat_occupancy='throat.occupancy_wp',
                                     pore_occupancy='pore.occupancy_wp',
                                     pdiffusive_conductance='pore.diffusive_conductance_wp',
                                     tdiffusive_conductance='throat.diffusive_conductance_wp')
            phys_moisture[key].models.add(propname='throat.conduit_conductance_dp',
                                     model=pm.multiphase.mixed_conductance,
                                     throat_occupancy='throat.occupancy_dp',
                                     pore_occupancy='pore.occupancy_dp',
                                     pdiffusive_conductance='pore.diffusive_conductance_dp',
                                     tdiffusive_conductance='throat.diffusive_conductance_dp')

            phys_moisture[key].models.regenerate()

            bounds = [['inlet', 'outlet']]
            pc1_wp = Pc_step_wp + dPc
            pc2_wp = Pc_step_wp - dPc
            pc1_dp = Pc_step_dp + dPc
            pc2_dp = Pc_step_dp - dPc

            for bound_increment in range(len(bounds)):
                alg_flow_wp[key] = pab.MoistureFlow(name='alg_flow_wp'+names[key], network=pn[key],
                                               phase=moisture[key])
                alg_flow_dp[key] = pab.MoistureFlow(name='alg_flow_dp'+names[key], network=pn[key],
                                               phase=moisture[key])

                BC1_pores = pn[key].pores(labels=bounds[bound_increment][0])
                BC2_pores = pn[key].pores(labels=bounds[bound_increment][1])

                # BC1
                alg_flow_wp[key].set_boundary_conditions(bctype='Dirichlet',
                                                    bcvalue=pc1_wp,
                                                    pores=BC1_pores)
                alg_flow_dp[key].set_boundary_conditions(bctype='Dirichlet',
                                                    bcvalue=pc1_dp,
                                                    pores=BC1_pores)

                # BC2
                alg_flow_wp[key].set_boundary_conditions(bctype='Dirichlet',
                                                    bcvalue=pc2_wp,
                                                    pores=BC2_pores)
                alg_flow_dp[key].set_boundary_conditions(bctype='Dirichlet',
                                                    bcvalue=pc2_dp,
                                                    pores=BC2_pores)

                # run algorithms with proper conduit conductance
                alg_flow_wp[key].run(conductance='conduit_conductance_wp',
                                quantity='pressure_wp')
                alg_flow_dp[key].run(conductance='conduit_conductance_dp',
                                quantity='pressure_dp')

                # calc effective permeabilities [s]
                eff_permeability_moisture_wp = alg_flow_wp[key].calc_eff_permeability(
                    conductance=phys_moisture[key]['throat.conduit_conductance_wp'])
                eff_permeability_moisture_dp = alg_flow_dp[key].calc_eff_permeability(
                    conductance=phys_moisture[key]['throat.conduit_conductance_dp'])

                # append permeability & flow values to the lists
                eff_perm_moisture_wp[str(bound_increment)].append(
                    eff_permeability_moisture_wp)
                eff_perm_moisture_dp[str(bound_increment)].append(
                    eff_permeability_moisture_dp)

                ctrl.purge_object(alg_flow_wp[key])
                ctrl.purge_object(alg_flow_dp[key])

                #% Verification: compare water occupancy
                # --------------------------------------
                Pc_p = alg_flow_wp[key]['pore.'+moisture[key]._name+'_pressure_wp']     # Pc result
                connected_pores = pn[key]['throat.conns']
                Pc_connected_pore = [[Pc_p[pair]] for pair in connected_pores]
                Pc_connected_pore = sp.array(Pc_connected_pore).reshape(
                                            (sp.shape(connected_pores)))
                Pc_t_result = sp.amin(Pc_connected_pore, axis=1)
    #            Pc_t_result = pn[key].interpolate_data(data=Pc_p)  # mean of 2 pores

                Tinvaded1 = water[key]['throat.occupancy_wp']
                Tinvaded2 = sp.float64(alg_wp[key]._t_cap<=Pc_step_wp)
                Tinvaded3 = sp.float64(alg_wp[key]._t_cap<=Pc_t_result)
                diff12 = sp.where(Tinvaded1!=Tinvaded2)
                diff23 = sp.where(Tinvaded3!=Tinvaded2)
                if sp.size(diff12):
                    print('Different12 water distribution at: ',
                          diff12)
                    print('Pc throat: ', alg_wp[key]._t_cap[diff12])
                    print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                    print('Pc step throat: ', Pc_t_result[diff12])
                    print('Pc step conn pores: ', Pc_connected_pore[diff12])

                if sp.size(diff23):
                    print('Different23 water distribution at: ',
                          diff23)
                    print('Pc throat: ', alg_wp[key]._t_cap[diff23])
                    print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                    print('Pc step throat: ', Pc_t_result[diff23])
                    print('Pc step conn pores: ', Pc_connected_pore[diff23])
                Ax = alg_flow_wp[key].A.dot(alg_flow_wp[key].X)
                b = alg_flow_wp[key].b.reshape(sp.shape(Ax))
                res = Ax - b
                norm_res = sp.linalg.norm(res, ord=2)
    #            rel_error = sp.linalg.norm(Ax)/sp.linalg.norm(b)
    #            print('Relative error: ', rel_error)
                print('Residual 2-norm: ', norm_res, '\n')
                max_norm_res = sp.amax(max_norm_res, norm_res)

        alg_flow_wp[key].store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                                 sat_moist=sat_wp_moist,
                                 w_sat=w_sat, k=eff_perm_moisture_wp['0'])
        alg_flow_dp[key].store_result(Pc=Pc_dp, sat=sat_dp, sat_surf=sat_dp_surf,
                                 sat_moist=sat_dp_moist,
                                 w_sat=w_sat, k=eff_perm_moisture_dp['0'])

        alg_flow_wp[key].calc_eff_permeability()
        alg_flow_wp[key].calc_abs_permeability()
        alg_flow_wp[key].calc_mD_permeability()
        alg_flow_dp[key].calc_eff_permeability()
        alg_flow_dp[key].calc_abs_permeability()
        alg_flow_dp[key].calc_mD_permeability()

        if plot:
            bpnm.Postprocessing.Plots.hysteresis(alg_flow_wp[key], alg_flow_dp[key])
    return (alg_flow_wp, alg_flow_dp)


