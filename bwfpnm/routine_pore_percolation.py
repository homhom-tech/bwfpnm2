# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:36:38 2016

@author: islah

This is a module which consists of routine functions for moisture storage and transfer estimations for topological network model

--> use throat.porelengths in conductivity calculation

Customised for Percolation class

"""
import bwfpnm as bpnm
import scipy as sp
from numpy.linalg import cond

ctrl = bpnm.Base.Controller()
ctrl.loglevel = 40

#%% Network
def network(NetClass=None, netname='net', dat_file=None, trimming=False,
            **kwargs):

    if dat_file is not None:
        netinput, geoinput, geomacro =\
            bpnm.Utilities.IO.load_data(dat_file)

    #% Create network object
    if NetClass is None:
        NetClass = bpnm.Network.Topology

    try:
        pn = NetClass(name=netname, **netinput)
    except:
        pn = NetClass(name=netname, **kwargs)

    #% Trim 'isolated' cluster
    if trimming:
        health = pn.check_network_health()
        trim_pores = health['trim_pores']
        # modify geoinput based on trimpore
        pn.trim_geom_data(geoinput, trim_pores)
        pn.trim(pores=trim_pores)

    try:
        return (pn, geoinput, geomacro)
    except:
        return pn


#%% Geometry
def geometry(NetInstance, GeoClass=None, geoname='geo',
             pores='all', throats='all', geoinput=None, **kwargs):

    ps = NetInstance.pores(pores)
    ts = NetInstance.throats(throats)

    if GeoClass is None:
        try:
            GeoClass = bpnm.Geometry.Topology
        except:
            GeoClass = bpnm.Geometry.Test_PoreThroat

    try:
        geo = GeoClass(network=NetInstance, pores=ps, throats=ts,
                       name=geoname, **geoinput)
    except:
        geo = GeoClass(network=NetInstance, pores=ps, throats=ts,
                       name=geoname, **kwargs)

    return geo


#%% Create phases: water + vapour = moisture
def phase(NetInstance, WaterClass=None, VapourClass=None,
          watername='water', vapourname='vapour', temperature=293.15,
          **kwargs):

    if WaterClass is None:
        WaterClass = bpnm.Phases.Water
    if VapourClass is None:
        VapourClass = bpnm.Phases.Vapour

    water = WaterClass(name=watername, network=NetInstance, **kwargs)
    vapour = VapourClass(name=vapourname, network=NetInstance, **kwargs)
    vapour['pore.water_density'] = water['pore.density']

    moisture = bpnm.Phases.GenericPhase(name='moisture', network=NetInstance,
                                        **kwargs)
    moisture['pore.temperature'] = temperature
    return (water, vapour, moisture)


#%% Create physics
def physic(NetInstance, WaterInstance, VapourInstance, MoistInstance,
           pores='all', throats='all',
           PhysWaterClass=None, PhysVapourClass=None, PhysMoistClass=None,
           physwatername='physics_water', physvapourname='physics_vapour',
           physmoistname='physics_moisture', param='shapefactor', **kwargs):

    if PhysWaterClass is None:
        PhysWaterClass = bpnm.Physics.Standard_Topology_pore
    if PhysVapourClass is None:
        PhysVapourClass = bpnm.Physics.GenericPhysics
    if PhysMoistClass is None:
        PhysMoistClass = bpnm.Physics.GenericPhysics

    ps = NetInstance.pores(pores)
    ts = NetInstance.throats(throats)
    # The physics for water is calculated beforehand since it is independent of
    # the ambient condition (rh or Pc)
    phys_water = PhysWaterClass(network=NetInstance, phase=WaterInstance,
                                pores=ps, throats=ts, param=param,
                                dynamic_data=True, name=physwatername)
    # The physics of vapour (and hence moisture) is calculated on the fly
    # since it is dependent to the ambient conditin (rh or Pc)
    phys_vapour = PhysVapourClass(network=NetInstance, phase=VapourInstance,
                                  pores=ps, throats=ts, dynamic_data=True,
                                  name=physvapourname)
    phys_moisture = PhysMoistClass(network=NetInstance, phase=MoistInstance,
                                   pores=ps, throats=ts, dynamic_data=True,
                                   name=physmoistname)

    return (phys_water, phys_vapour, phys_moisture)


#%% Create algorithm: Wetting & Drying Percolation -> only use element's Pc_wp
def algorithm(NetIns, WaterIns, VapourIns, Class=None, case='wetting',
              npts=None, dp=False, trapping=False,
              inlet_imb=None, inlet_dry=None, outlet=None, **kwargs):

    if Class is None:
        Class = bpnm.Algorithms.Percolation
    alg_wp = Class(network=NetIns, phase_wet=WaterIns, phase_dry=VapourIns,
                   name='percolation', **kwargs)
    if inlet_dry is None:
        inlet_dry = NetIns['pore.inlet'] + NetIns['pore.outlet']    # bool
    if inlet_imb is None:
        inlet_imb = NetIns['pore.inlet']    # bool

    if case=='wetting':
        alg_wp.run_wetting(prop_name='wetting_inv_seq')
        if trapping:
            alg_wp.evaluate_trapping_wetting(p_outlets=outlet, mode='clone')
            alg_wp.copy_results(prop_names=[case+'_trapped_pc'])
    elif case=='imbibition':
        alg_wp.set_inlets_imbibition(pores=inlet_imb)
        alg_wp.run_imbibition(prop_name='imbibition_inv_seq')
        if trapping:
            alg_wp.evaluate_trapping_imbibition(p_outlets=outlet, mode='clone')

    # copy result to wet phase
    alg_wp.copy_results(prop_names=[case+'_inv_seq', case+'_inv_pc',
                              case+'_inv_sat'])
    if trapping:
        alg_wp.copy_results(prop_names=[case+'_inv_seq_trapping',
                                        case+'_inv_pc_trapping',
                                        case+'_inv_sat_trapping',
                                        case+'_trapped_seq'])
    # if drying from wetting/imbibition
    if dp:
        alg_dp = Class(network=NetIns, phase_wet=WaterIns, phase_dry=VapourIns,
                       name='percolation_dp', **kwargs)
        if trapping:
            pinv = alg_wp['pore.'+case+'_inv_pc_trapping']
            tinv = alg_wp['throat.'+case+'_inv_pc_trapping']
        else:
            pinv = alg_wp['pore.'+case+'_inv_pc']
            tinv = alg_wp['throat.'+case+'_inv_pc']

        alg_dp.setup_drying(p_inv_pc=pinv, t_inv_pc=tinv)
        alg_dp.run_drying(inv_site=inlet_dry)
        # copy result to wet phase
        case = 'drying'
        alg_dp.copy_results(prop_names=[case+'_inv_seq', case+'_inv_pc',
                                        case+'_inv_sat'])
    else:
        alg_dp = 0

    return (alg_wp, alg_dp)


#%% create moisture content (w) array
def moist_content(geo, alg_wp, water_density, case='wetting',
                  trapping=False, alg_dp=None,
                  v_mat=None, porosity=None, **kwargs):

    v_pore = sp.sum(geo['pore.volume']) + sp.sum(geo['throat.volume'])

    if v_mat is None:
        v_mat = geo._net._macro_Lx * geo._net._macro_Ly * geo._net._macro_Lz

    porosity = v_pore/v_mat

    try:
        w_sat = porosity*water_density
    except:
        print('error: either volume of bulk material of porosity is required!')
        return

    alg_wp['pore.'+case+'_inv_w'] = alg_wp['pore.'+case+'_inv_sat']*w_sat
    alg_wp['throat.'+case+'_inv_w'] = alg_wp['throat.'+case+'_inv_sat']*w_sat
    if trapping:
        prop = case+'_inv_w_trapping'
        prop_sat = case+'_inv_sat_trapping'
        alg_wp['pore.'+prop] = alg_wp['pore.'+prop_sat]*w_sat
        alg_wp['throat.'+prop] = alg_wp['throat.'+prop_sat]*w_sat

        alg_wp.copy_results(prop_names=[case+'_inv_w_trapping'])
    alg_wp.copy_results(prop_names=[case+'_inv_w'])

    if alg_dp is not None:
        case = 'drying'
        prop = case+'_inv_w'
        alg_dp['pore.'+prop] = alg_dp['pore.'+case+'_inv_sat']*w_sat
        alg_dp['throat.'+prop] = alg_dp['throat.'+case+'_inv_sat']*w_sat

        alg_dp.copy_results(prop_names=[case+'_inv_w'])

    return (alg_wp, alg_dp, w_sat, porosity)


#%% Calculate permeability for each Pc value
def permeability(pn, alg_wp, alg_dp, water, vapour, moisture,
                 phys_vapour, phys_moisture, w_sat, case='wetting',
                 trapping=False,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc =1, diffusion=True,
                 **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat_wp, sat_dp = [], []
    sat_wp_surf, sat_dp_surf = [], []
    sat_wp_moist, sat_dp_moist = [], []
#    max_norm_res = 0
    eff_perm_moisture_wp = {'0': [], '1': [], '2': []}
    eff_perm_moisture_dp = {'0': [], '1': [], '2': []}

    p_volumes = pn['pore.volume']
    t_volumes = pn['throat.volume']
    volume_total = sum(p_volumes) + sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
#    Pc_wp = sp.around(Pc_wp, 3)
    Pc_dp = Pc_wp[::-1]

    occupancy = 'occupancy_'+case
    if trapping:
        occupy = occupancy + '_trapping'
    else:
        occupy = occupancy

    for Pc_step_wp, Pc_step_dp in list(zip(Pc_wp, Pc_wp)):
        alg_wp.return_results(Pc=Pc_step_wp, cases=[case],
                              occupancies=[occupancy], trapping=trapping)
        alg_dp.return_results(Pc=Pc_step_dp, cases=['drying'],
                              occupancies=['occupancy_drying'],
                              trapping=False)
        if (Pc_step_dp+120094) <50000:
            pass

        p_occ_wp = water['pore.'+occupy]
        t_occ_wp = water['throat.'+occupy]
        p_occ_dp = water['pore.occupancy_drying']
        t_occ_dp = water['throat.occupancy_drying']

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
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_wp,
                                   pore_occupancy='pore.'+occupy)
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_wp,
                                   throat_occupancy='throat.'+occupy)
            phys_vapour.models.add(propname='pore.surface_thickness_drying',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_drying')
            phys_vapour.models.add(propname='throat.surface_thickness_drying',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_drying')

            ## Inner circular tube
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.'+occupy,
                                   film_thickness='pore.surface_thickness_wp')
            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.'+occupy,
                                   film_thickness='throat.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_drying',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.occupancy_drying',
                                   film_thickness='pore.surface_thickness_drying')
            phys_vapour.models.add(propname='throat.surface_volume_drying',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.occupancy_drying',
                                   film_thickness='throat.surface_thickness_drying')


            volume_p_wp = sum(phys_vapour['pore.surface_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])
            volume_p_dp = sum(phys_vapour['pore.surface_volume_drying'])
            volume_t_dp = sum(phys_vapour['throat.surface_volume_drying'])

            sat_surf_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_surf_dp = (volume_p_dp + volume_t_dp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)
            sat_dp_surf.append(sat_surf_dp)

            sat_wp[-1] += sat_surf_wp   # update total saturations
            sat_dp[-1] += sat_surf_dp

            print('pthickness wp: ',
                  phys_vapour['pore.surface_thickness_wp'][0])
            print('pradius wp: ',
                  pn['pore.diameter'][0]/2)
            print('pvol surf wp: ',
                  phys_vapour['pore.surface_volume_wp'][0])
            print('pvol wp: ',
                  p_volumes[0])
            print('psat wp: ',
                  phys_vapour['pore.surface_volume_wp'][0]/p_volumes[0])

        if moist_volume:
            phys_vapour.models.add(propname='pore.moist_volume_wp',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='throat.moist_volume_wp',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.moist_volume_drying',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_drying')
            phys_vapour.models.add(propname='throat.moist_volume_drying',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_drying')


            volume_p_wp = sum(phys_vapour['pore.moist_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.moist_volume_wp'])
            volume_p_dp = sum(phys_vapour['pore.moist_volume_drying'])
            volume_t_dp = sum(phys_vapour['throat.moist_volume_drying'])

            sat_moist_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_moist_dp = (volume_p_dp + volume_t_dp)/(volume_total)
            sat_wp_moist.append(sat_moist_wp)
            sat_dp_moist.append(sat_moist_dp)

            sat_wp[-1] += sat_moist_wp   # update total saturations
            sat_dp[-1] += sat_moist_dp

            print('moist vol: ', phys_vapour['throat.moist_volume_wp'])
            print('moist sat: ',
                  phys_vapour['throat.moist_volume_wp']/t_volumes)

        # Update vapour permeability for ALL pores & throats for wp & dp
        # g = f(Pv(rh(pc)))
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen,
                               film_thickness='throat.surface_thickness_wp')
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_wp, knudsen=knudsen,
                                film_thickness='pore.surface_thickness_wp')

        phys_vapour.models.add(propname='throat.diffusive_conductance_drying',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_dp, knudsen=knudsen,
                               film_thickness='throat.surface_thickness_drying')
        phys_vapour.models.add(propname='throat.diffusive_conductance_drying_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_dp, knudsen=knudsen,
                                film_thickness='pore.surface_thickness_drying')
        phys_vapour.models.regenerate()

        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.'+occupy,
                                 pore_occupancy='pore.'+occupy,
                                 pdiffusive_conductance='throat.diffusive_conductance_wp_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp',
                                 diffusion=diffusion)
        phys_moisture.models.add(propname='throat.conduit_conductance_drying',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.occupancy_drying',
                                 pore_occupancy='pore.occupancy_drying',
                                 pdiffusive_conductance='throat.diffusive_conductance_drying_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_drying',
                                 diffusion=diffusion)

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp + dPc
        pc2_wp = Pc_step_wp - dPc
        pc1_dp = Pc_step_dp + dPc
        pc2_dp = Pc_step_dp - dPc

        for bound_increment in range(len(bounds)):
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            alg_flow_dp = pab.MoistureFlow(name='alg_flow_dp', network=pn,
                                           phase=moisture)

            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # BC1
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_wp,
                                                pores=BC1_pores)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_dp,
                                                pores=BC1_pores)

            # BC2
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_wp,
                                                pores=BC2_pores)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_dp,
                                                pores=BC2_pores)

            # run algorithms with proper conduit conductance
            alg_flow_wp.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')
            alg_flow_dp.run(conductance='conduit_conductance_drying',
                                quantity='pressure_drying')
#            try:
#                alg_flow_dp.run(conductance='conduit_conductance_drying',
#                                quantity='pressure_drying',
#                                iterative_solver='cg',
#                                x0=x0) #, tol=1e-10, maxiter=1000,xtype, M,callback)
#                x0 = alg_flow_dp.X
#                print('iterative solver has been run at Pc = ', Pc_step_dp)
#            except:
#                alg_flow_dp.run(conductance='conduit_conductance_drying',
#                                quantity='pressure_drying')
#                x0 = alg_flow_dp.X

            # calc effective permeabilities [s]
            eff_permeability_moisture_wp = alg_flow_wp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            eff_permeability_moisture_dp = alg_flow_dp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_drying'])

            # append permeability & flow values to the lists
            eff_perm_moisture_wp[str(bound_increment)].append(
                eff_permeability_moisture_wp)
            eff_perm_moisture_dp[str(bound_increment)].append(
                eff_permeability_moisture_dp)

            ctrl.purge_object(alg_flow_wp)
            ctrl.purge_object(alg_flow_dp)

#            #% Verification: compare water occupancy
#            # --------------------------------------
#            Pc_p = alg_flow_wp['pore.moisture_pressure_'+case]     # Pc result
#            connected_pores = pn['throat.conns']
#            Pc_connected_pore = [[Pc_p[pair]] for pair in connected_pores]
#            Pc_connected_pore = sp.array(Pc_connected_pore).reshape(
#                                        (sp.shape(connected_pores)))
#            Pc_t_result = sp.amin(Pc_connected_pore, axis=1)
##            Pc_t_result = pn.interpolate_data(data=Pc_p)  # mean of 2 pores
#
#            Tinvaded1 = water['throat.occupancy_'+case]
#            Tinvaded2 = sp.float64(alg_wp._t_cap<=Pc_step_wp)
#            Tinvaded3 = sp.float64(alg_wp._t_cap<=Pc_t_result)
#            diff12 = sp.where(Tinvaded1!=Tinvaded2)
#            diff23 = sp.where(Tinvaded3!=Tinvaded2)
#            if sp.size(diff12):
#                print('Different12 water distribution at: ',
#                      diff12)
#                print('Pc throat: ', alg_wp._t_cap[diff12])
#                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
#                print('Pc step throat: ', Pc_t_result[diff12])
#                print('Pc step conn pores: ', Pc_connected_pore[diff12])
#
#            if sp.size(diff23):
#                print('Different23 water distribution at: ',
#                      diff23)
#                print('Pc throat: ', alg_wp._t_cap[diff23])
#                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
#                print('Pc step throat: ', Pc_t_result[diff23])
#                print('Pc step conn pores: ', Pc_connected_pore[diff23])
#            Ax = alg_flow_wp.A.dot(alg_flow_wp.X)
#            b = alg_flow_wp.b.reshape(sp.shape(Ax))
#            res = Ax - b
##            norm_res = sp.linalg.norm(res, ord=2)
##            print('Residual 2-norm: ', norm_res, '\n')
##            max_norm_res = sp.amax(max_norm_res, norm_res)

    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=eff_perm_moisture_wp['0'])
    alg_flow_dp.store_result(Pc=Pc_wp, sat=sat_dp, sat_surf=sat_dp_surf,
                             sat_moist=sat_dp_moist,
                             w_sat=w_sat, k=eff_perm_moisture_dp['0'])

    alg_flow_wp.calc_abs_permeability()
    alg_flow_wp.calc_mD_permeability()
    alg_flow_dp.calc_abs_permeability()
    alg_flow_dp.calc_mD_permeability()

    if plot:
        bpnm.Postprocessing.Plots.hysteresis(alg_flow_wp, alg_flow_dp)
    return (alg_flow_wp, alg_flow_dp)


def permeability2(pn, alg_wp, alg_dp, water, vapour, moisture,
                 phys_vapour, phys_moisture, w_sat, case='wetting',
                 trapping=False,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc =1, diffusion=True,
                 refine=False, debug=False,
                 **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    pcs = []
    sat_wp, sat_dp = [], []
    sat_wp_surf, sat_dp_surf = [], []
    sat_wp_moist, sat_dp_moist = [], []
    eff_perm_moisture_wp = {'0': [], '1': [], '2': []}
    eff_perm_moisture_dp = {'0': [], '1': [], '2': []}

    p_volumes = pn['pore.volume']
    t_volumes = pn['throat.volume']
    volume_total = sum(p_volumes) + sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
#    Pc_wp = sp.around(Pc_wp, 3)
    Pc_list = sp.copy(Pc_wp)
    Pc_wp = list(Pc_wp[::-1])
#    Pc_dp = Pc_wp[::-1]

    occupancy = 'occupancy_'+case
    if trapping:
        occupy = occupancy + '_trapping'
    else:
        occupy = occupancy
    i = 0
    while len(Pc_wp) > 0:
        Pc_step_wp = Pc_wp.pop()
        Pc_step_dp = Pc_step_wp

#        try:
#            Pc_step_dp = Pc_dp.pop()
#        except:
#            Pc_step_dp = Pc_step_wp
        if i==14:
            pass
        alg_wp.return_results(Pc=Pc_step_wp, cases=[case],
                              occupancies=[occupancy], trapping=trapping)
        alg_dp.return_results(Pc=Pc_step_dp, cases=['drying'],
                              occupancies=['occupancy_drying'],
                              trapping=False)

        p_occ_wp = water['pore.'+occupy]
        t_occ_wp = water['throat.'+occupy]
        p_occ_dp = water['pore.occupancy_drying']
        t_occ_dp = water['throat.occupancy_drying']

        check_occu_wp = _check_occupancy(p_occ_wp, t_occ_wp, 'wetting')
        check_occu_dp = _check_occupancy(p_occ_dp, t_occ_dp, 'drying')

#        if i in [14, 15]:
#            p_occ_wp = water['pore.occupancy_drying']
#            t_occ_wp = water['throat.occupancy_drying']

        volume_p_wp = sum(p_occ_wp*p_volumes)
        volume_t_wp = sum(t_occ_wp*t_volumes)
        volume_p_dp = sum(p_occ_dp*p_volumes)
        volume_t_dp = sum(t_occ_dp*t_volumes)

        saturation_wp = (volume_p_wp + volume_t_wp)/(volume_total)
        saturation_dp = (volume_p_dp + volume_t_dp)/(volume_total)

        sat_wp.append(saturation_wp)
        sat_dp.append(saturation_dp)
        pcs.append(Pc_step_wp)

        # Surface adsorption: p/t.surface_thickness & surface_volume
        if surface_ad:
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_wp,
                                   pore_occupancy='pore.'+occupy)
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_wp,
                                   throat_occupancy='throat.'+occupy)
            phys_vapour.models.add(propname='pore.surface_thickness_drying',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_drying')
            phys_vapour.models.add(propname='throat.surface_thickness_drying',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_drying')

            ## Inner circular tube
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.'+occupy,
                                   film_thickness='pore.surface_thickness_wp')
            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.'+occupy,
                                   film_thickness='throat.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_drying',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.occupancy_drying',
                                   film_thickness='pore.surface_thickness_drying')
            phys_vapour.models.add(propname='throat.surface_volume_drying',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.occupancy_drying',
                                   film_thickness='throat.surface_thickness_drying')


            volume_p_wp = sum(phys_vapour['pore.surface_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])
            volume_p_dp = sum(phys_vapour['pore.surface_volume_drying'])
            volume_t_dp = sum(phys_vapour['throat.surface_volume_drying'])

            sat_surf_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_surf_dp = (volume_p_dp + volume_t_dp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)
            sat_dp_surf.append(sat_surf_dp)

            sat_wp[-1] += sat_surf_wp   # update total saturations
            sat_dp[-1] += sat_surf_dp

        # Update vapour permeability for ALL pores & throats for wp & dp
        # g = f(Pv(rh(pc)))
        # for the same pcwp & pcdp, diffusive_conductance of wp & dp are the same
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen,
                               film_thickness='throat.surface_thickness_wp')
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_wp, knudsen=knudsen,
                                film_thickness='pore.surface_thickness_wp')

        phys_vapour.models.add(propname='throat.diffusive_conductance_drying',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_dp, knudsen=knudsen,
                               film_thickness='throat.surface_thickness_drying')
        phys_vapour.models.add(propname='throat.diffusive_conductance_drying_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_dp, knudsen=knudsen,
                                film_thickness='pore.surface_thickness_drying')
        phys_vapour.models.regenerate()

        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.'+occupy,
                                 pore_occupancy='pore.'+occupy,
                                 pdiffusive_conductance='throat.diffusive_conductance_wp_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp',
                                 diffusion=diffusion)
        phys_moisture.models.add(propname='throat.conduit_conductance_drying',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.occupancy_drying',
                                 pore_occupancy='pore.occupancy_drying',
                                 pdiffusive_conductance='throat.diffusive_conductance_drying_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_drying',
                                 diffusion=diffusion)
#        if i in [14, 15]:
#            phys_moisture.models.add(propname='throat.conduit_conductance_wp',
#                                 model=pm.multiphase.mixed_conductance_pore,
#                                 throat_occupancy='throat.occupancy_drying',
#                                 pore_occupancy='pore.occupancy_drying',
#                                 pdiffusive_conductance='throat.diffusive_conductance_drying_pore',
#                                 tdiffusive_conductance='throat.diffusive_conductance_drying',
#                                 diffusion=diffusion)

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp + dPc
        pc2_wp = Pc_step_wp - dPc
        pc1_dp = Pc_step_dp + dPc
        pc2_dp = Pc_step_dp - dPc

        kmin = phys_moisture['throat.conduit_conductance_drying'].min()
        kmax = phys_moisture['throat.conduit_conductance_drying'].max()
        print('Pc: {}'.format(Pc_step_wp))
#        print('Pc: {}, kmin: {}, kmax: {}'.format(Pc_step_wp, kmin, kmax))

        for bound_increment in range(len(bounds)):
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            alg_flow_dp = pab.MoistureFlow(name='alg_flow_dp', network=pn,
                                           phase=moisture)

            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # BC1
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_wp,
                                                pores=BC1_pores)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_dp,
                                                pores=BC1_pores)

            # BC2
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_wp,
                                                pores=BC2_pores)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_dp,
                                                pores=BC2_pores)

            # run algorithms with proper conduit conductance
            alg_flow_wp.run2(conductance='conduit_conductance_wp',
                            quantity='pressure_wp', **kwargs)
            alg_flow_dp.run2(conductance='conduit_conductance_drying',
                                quantity='pressure_drying', **kwargs)

            # calc effective permeabilities [s]
            eff_permeability_moisture_wp = alg_flow_wp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            eff_permeability_moisture_dp = alg_flow_dp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_drying'])

            # append permeability & flow values to the lists
            eff_perm_moisture_wp[str(bound_increment)].append(
                eff_permeability_moisture_wp)
            eff_perm_moisture_dp[str(bound_increment)].append(
                eff_permeability_moisture_dp)

            alg_flow_wp.return_rate(case='_wetting')
            alg_flow_dp.return_rate(case='_drying')

            if debug:
                rate1 = moisture['throat.rate_drying']
                deltaP1 = moisture['throat.delta_pressure_drying']
                gm1 = alg_flow_dp['throat.conductance']
                X1 = alg_flow_dp.X
                A1 = alg_flow_dp.A
                b1 = alg_flow_dp.b
                bb = A1*X1.reshape(b1.shape)
                err = sp.linalg.norm(bb-b1)
                berr = sp.amax(sp.absolute(bb[b1==0]))
                try:
                    tlocrate = sp.where(rate1<rate0)[0]
                    tlocg = sp.where(gm1<gm0)[0]
                    tlocdP = sp.where(sp.absolute(deltaP1)<sp.absolute(deltaP0))[0]
                    gdiff = gm1[tlocrate] < gm0[tlocrate]
                    dPdiff = deltaP1[tlocrate] < deltaP0[tlocrate]
                    if tlocrate:
                        pass
                except:
                    pass
                rate0, gm0, deltaP0 = rate1, gm1, deltaP1
                X0, A0, b0 = X1, A1, b1

            ctrl.purge_object(alg_flow_wp)
            ctrl.purge_object(alg_flow_dp)
        # Grid refinement
        if len(Pc_wp)==0 and refine:
            arg = sp.argsort(pcs)
            pcs = sp.array(pcs)[arg]
            sat_wp = sp.array(sat_wp)[arg]
            eff_perm_moisture_wp['0'] = sp.array(eff_perm_moisture_wp['0'])[arg]
            sat_dp = sp.array(sat_dp)[arg]
            eff_perm_moisture_dp['0'] = sp.array(eff_perm_moisture_dp['0'])[arg]
            Pc_list = Pc_list[arg]

            if surface_ad:
                sat_wp_surf = sp.array(sat_wp_surf)[arg]
                sat_dp_surf = sp.array(sat_dp_surf)[arg]

            dw = sp.absolute(sp.diff(sat_wp))
            dw2 = sp.absolute(sp.diff(sat_dp))
            tag = sp.where(dw>0.2)[0]
            tag2 = sp.where(dw2>0.2)[0]
            tag = sp.unique(sp.r_[tag, tag2])
            tagged = sp.array([Pc_list[tag], Pc_list[tag+1]])
            fine = -sp.power(10, sp.mean(sp.log10(-tagged), axis=0))
            if len(fine)>0:
                Pc_wp.extend(fine)
                Pc_list = sp.concatenate((Pc_list,fine[::-1]))

                pcs = list(pcs)
                sat_wp = list(sat_wp)
                sat_wp_surf = list(sat_wp_surf)
                eff_perm_moisture_wp['0'] = list(eff_perm_moisture_wp['0'])
                sat_dp = list(sat_dp)
                sat_dp_surf = list(sat_dp_surf)
                eff_perm_moisture_dp['0'] = list(eff_perm_moisture_dp['0'])
        i += 1
        A = sp.absolute(alg_flow_dp.A.data)
        print('Amin: {}, Amax: {}'.format(A.min(), A.max()))
#        Ad = alg_flow_dp.A
#        bd = alg_flow_dp.b
#        xd = alg_flow_dp.X
#        xi, istop, itn, normr, normar, normA, condA, normx = sp.sparse.linalg.lsmr(Ad, bd, atol=0, btol=0, conlim=0, maxiter=1000)
#        xii, istop2, itn2, r1norm, r2norm, anorm, acond, arnorm, xnorm, var = sp.sparse.linalg.lsqr(Ad, bd, damp=0.0, atol=0, btol=0, conlim=0, iter_lim=1000, show=False, calc_var=False)
#        xerr = sp.linalg.norm(xd-xi, 2)
#        Adense = Ad.todense()
#        conda = cond(Adense)
        # condition number = sigma.max/sigma.min, SVD
#        U, s, Vh = sp.sparse.linalg.svds(Ad, k=Ad.shape[0]-1, maxiter=None)
#        U, s, Vh = sp.linalg.svd(Ad, full_matrices=False)
#        print('conda: {}'.format(conda))
#        print('Cond_lsq: {}, ||xd-xi||: {}\n'.format(condA, xerr))
#
#        if conda >= 1./sp.finfo(sp.float64).eps:
#            print('WARNING!! Ill-conditioned system: huge condition number\n')

    arg = sp.argsort(pcs)
    pcwp = sp.array(pcs)[arg]
    sat_wp = sp.array(sat_wp)[arg]
    k = sp.array(eff_perm_moisture_wp['0'])[arg]
    if surface_ad:
        sat_wp_surf = sp.array(sat_wp_surf)[arg]
    alg_flow_wp.store_result(Pc=pcwp, sat=sat_wp,
                             sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=k)

    pcs = sp.array(pcs)[arg]
    sat_dp = sp.array(sat_dp)[arg]
    k = sp.array(eff_perm_moisture_dp['0'])[arg]
    if surface_ad:
        sat_dp_surf = sp.array(sat_dp_surf)[arg]
    alg_flow_dp.store_result(Pc=pcs, sat=sat_dp,
                             sat_surf=sat_dp_surf,
                             sat_moist=sat_dp_moist,
                             w_sat=w_sat, k=k)

    alg_flow_wp.calc_abs_permeability()
    alg_flow_wp.calc_mD_permeability()
    alg_flow_dp.calc_abs_permeability()
    alg_flow_dp.calc_mD_permeability()

    if plot:
        bpnm.Postprocessing.Plots.hysteresis(alg_flow_wp, alg_flow_dp)
    return (alg_flow_wp, alg_flow_dp)

def _check_occupancy(p_occ, t_occ, case, conns):
    if case=='wetting':
        pass
    elif case=='drying':
        t_pocc = p_occ[conns]
        psum = sp.sum(t_pocc, axis=0)

# %% A group of steps to calculate moisture retention curve
def moisture_storage(plot=True, **kwargs):
    pn, geoinput, geomacro = network(
        dat_file='test/smallBenth.p',
        trimming=True)
    geo = geometry(pn, geoinput=geoinput)
    water, vapour, moisture = phase(pn)
    phys_water, phys_vapour, phys_moisture = physic(pn, water, vapour,
                                                    moisture)
    alg_wp, alg_dp = algorithm(pn, water, vapour)
    v_mat = sp.prod(geomacro[b'size'])
    alg_wp, alg_dp, w_sat, porosity = moist_content(geo, alg_wp, alg_dp,
                                                    water['pore.density'][0],
                                                    v_mat=v_mat)
    if plot:
        plot_moist_retention(alg_wp)

    return(pn, geo, water, vapour, moisture, phys_water, phys_vapour,
           phys_moisture, alg_wp, alg_dp, w_sat, porosity)


# %% Main
if __name__ == '__main__':
    (pn, geo, water, vapour, moisture, phys_water, phys_vapour, phys_moisture,
     alg_wp, alg_dp, w_sat, porosity) = moisture_storage(plot=False)
    alg_flow_wp, alg_flow_dp = permeability(pn, alg_wp, alg_dp,
                                            water, vapour, moisture,
                                            phys_vapour, phys_moisture,
                                            num_seq=10, w_sat=w_sat)
