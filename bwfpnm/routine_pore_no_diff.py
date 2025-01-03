# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:50:31 2015

@author: islah


This is a module which consists of routine functions for moisture storage and transfer estimations for topological network model

--> use throat.porelengths in conductivity calculation

"""
import bwfpnm as bpnm
import scipy as sp
import OpenPNM

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
          watername='water', vapourname='vapour', props=None,
          **kwargs):

    if WaterClass is None:
        WaterClass = bpnm.Phases.Water
    if VapourClass is None:
        VapourClass = bpnm.Phases.Vapour

    water = WaterClass(name=watername, network=NetInstance, props=props['water'],
                       **kwargs)
    vapour = VapourClass(name=vapourname, network=NetInstance, props=props['vapour'],
                         **kwargs)
    vapour['prop.water_density'] = water['prop.density']

    moisture = bpnm.Phases.GenericPhase(name='moisture', network=NetInstance,
                                        props=props['moisture'], **kwargs)
#    moisture['prop.temperature'] = temperature
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
                                  pores=ps, throats=ts,
                                  dynamic_data=True, name=physvapourname)
    phys_moisture = PhysMoistClass(network=NetInstance, phase=MoistInstance,
                                   pores=ps, throats=ts,
                                   dynamic_data=True, name=physmoistname)

    return (phys_water, phys_vapour, phys_moisture)


#%% Create algorithm: Wetting & Drying Percolation -> only use element's Pc_wp
def algorithm(NetIns, WaterIns, VapourIns, AdsorpClass=None, DesorpClass=None,
              npts=None, inv_points=None, dp=True, air_entrapment=False,
              inv_sites=None, **kwargs):

    if AdsorpClass is None:
        AdsorpClass = bpnm.Algorithms.WettingPercolation
    alg_wp = AdsorpClass(network=NetIns,
                         invading_phase=WaterIns, defending_phase=VapourIns,
                         name='WettingPercolation', **kwargs)
    inv_sites = NetIns['pore.inlet'] + NetIns['pore.outlet']    # bool
    alg_wp.run(inlets=inv_sites, npts=npts, inv_points=inv_points,
               access_limited=air_entrapment)

    if dp:
        if DesorpClass is None:
            DesorpClass = bpnm.Algorithms.DryingPercolation
        if npts is None:
            npts_dp = 10
        else:
            npts_dp = npts
        alg_dp = DesorpClass(network=NetIns,
                             invading_phase=VapourIns, defending_phase=WaterIns,
                             name='DryingPercolation', **kwargs)
        alg_dp.run(inlets=inv_sites, npts=npts_dp, inv_points=inv_points)
    else:
        alg_dp = 0

    return (alg_wp, alg_dp)


#%% create moisture content (w) array
def moist_content(geo, alg_wp, alg_dp, water_density,
                  v_mat=None, porosity=None, case='imbibition', **kwargs):

    try:
        v_pore = sp.sum(geo['pore.volume']) + sp.sum(geo['throat.volume'])
    except:
        v_pore = sp.sum(geo['pore.volume'])

    if v_mat is None:
        v_mat = geo._net._macro_Lx * geo._net._macro_Ly * geo._net._macro_Lz

    porosity = v_pore/v_mat

    try:
        w_sat = porosity*water_density
    except:
        print('error: either volume of bulk material of porosity is required!')
        return

    alg_wp['pore.'+case+'_inv_w'] = alg_wp['pore.'+case+'_inv_sat']*w_sat
    try:
        alg_dp['pore.'+case+'_inv_w'] = alg_dp['pore.'+case+'_inv_sat']*w_sat
    except:
        pass
    try:
        alg_wp['throat.'+case+'_inv_w'] = alg_wp['throat.'+case+'_inv_sat']*w_sat
        alg_dp['throat.'+case+'_inv_w'] = alg_dp['throat.'+case+'_inv_sat']*w_sat
    except:
        pass

    return (alg_wp, alg_dp, w_sat, porosity)


#%% Postprocessing: Plot
def plot_moist_retention(alg_wp, sat=None, Pc=None, **kwargs):

    if sat is not None:
        alg_wp.return_results(sat=sat)
    elif Pc is not None:
        alg_wp.return_results(Pc=Pc)
    else:
        alg_wp.return_results(sat=.5)
#     bpnm.Postprocessing.Plots.moisture_retention(alg_wp, w_sat)
    bpnm.Postprocessing.Plots.wetting_curves(alg_wp)


#%% Import for Paraview
def create_vtk(NetInstance, listPhases, filename='wetting.vtp', **kwargs):
     import bwfpnm.Utilities.IO as io
     io.VTK.save(network=NetInstance, filename=filename, phases=listPhases)
     # the new way in the new OpenPNM version
     # ctrl.save(), or
     # ctrl.save_object(alg_flow_wp, filename=)


#%% Calculate permeability for each Pc value
def permeability(pn, alg_wp, alg_dp, water, vapour, moisture,
                 phys_water, phys_vapour, phys_moisture, w_sat,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc =1, phase='water',
                 case='imbibition', **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat_wp, sat_dp = [], []
    sat_wp_surf, sat_dp_surf = [], []
    sat_wp_moist, sat_dp_moist = [], []
    span_air, span_water = [], []
#    max_norm_res = 0
    eff_perm_moisture_wp = {'0': [], '1': [], '2': []}
    eff_perm_moisture_dp = {'0': [], '1': [], '2': []}
    eff_perm_water_wp = {'0': [], '1': [], '2': []}
    eff_perm_vapour_wp = {'0': [], '1': [], '2': []}

    p_volumes = pn['pore.volume']
    t_volumes = pn['throat.volume']
    volume_total = sum(p_volumes) + sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
    Pc_wp = sp.around(Pc_wp, 3)
    Pc_dp = Pc_wp[::-1]

    for Pc_step_wp, Pc_step_dp in list(zip(Pc_wp, Pc_dp)):
        alg_wp.return_results(Pc=Pc_step_wp, occupancies=['occupancy_wp'], cases=[case])
        alg_dp.return_results(Pc=Pc_step_dp, occupancies=['occupancy_dp'], cases=['drying_'+case])

        p_occ_wp = water['pore.occupancy_wp']
        t_occ_wp = water['throat.occupancy_wp']
        p_occ_dp = water['pore.occupancy_dp']
        t_occ_dp = water['throat.occupancy_dp']

        #        if not t_occ_wp.min():
        tmask = ~sp.bool_(t_occ_wp)
        (aSpanBool, aSpanCl) = pn.span_existence(mask=tmask, return_isolated=False)
        if aSpanBool:
            print('Spanning air cluster exists at lPc: %.3f'
                    % sp.log10(-Pc_step_wp))
            span_air.append(Pc_step_wp)
        else:
            print('The biggest size of vapour cluster: %.3f' % sp.size(aSpanCl))

#        if p_occ_wp.max():
        (SpanBool, SpanCl) = pn.span_existence(mask=p_occ_wp, return_isolated=False)
        if SpanBool:
            print('Spanning water cluster exists at lPc: %.3f'
                   % sp.log10(-Pc_step_wp))
            span_water.append(Pc_step_wp)
        else:
            print('The biggest size of water cluster: %.3f'
                   % sp.size(SpanCl))

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
                                   model=pm.surface_adsorption.stat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.stat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.surface_thickness_dp',
                                   model=pm.surface_adsorption.stat_thickness,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_dp')
            phys_vapour.models.add(propname='throat.surface_thickness_dp',
                                   model=pm.surface_adsorption.stat_thickness,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_dp')

            ## Inner circular tube
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume)
            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume)
            phys_vapour.models.add(propname='pore.surface_volume_dp',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.occupancy_dp',
                                   film_thickness='pore.surface_thickness_dp')
            phys_vapour.models.add(propname='throat.surface_volume_dp',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.occupancy_dp',
                                   film_thickness='throat.surface_thickness_dp')


            volume_p_wp = sum(phys_vapour['pore.surface_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])
            volume_p_dp = sum(phys_vapour['pore.surface_volume_dp'])
            volume_t_dp = sum(phys_vapour['throat.surface_volume_dp'])

            sat_surf_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_surf_dp = (volume_p_dp + volume_t_dp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)
            sat_dp_surf.append(sat_surf_dp)

            sat_wp[-1] += sat_surf_wp   # update total saturations
            sat_dp[-1] += sat_surf_dp

        # Update vapour permeability for ALL pores & throats for wp & dp
        # g = f(Pv(rh(pc)))
#        if phase=='vapour':
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_wp, knudsen=knudsen)

        phys_vapour.models.add(propname='throat.diffusive_conductance_dp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_dp, knudsen=knudsen)
        phys_vapour.models.add(propname='throat.diffusive_conductance_dp_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_dp, knudsen=knudsen)
        phys_vapour.models.regenerate()

        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.occupancy_wp',
                                 pore_occupancy='pore.occupancy_wp',
                                 pdiffusive_conductance='throat.diffusive_conductance_wp_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp')
        phys_moisture.models.add(propname='throat.conduit_conductance_dp',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.occupancy_dp',
                                 pore_occupancy='pore.occupancy_dp',
                                 pdiffusive_conductance='throat.diffusive_conductance_dp_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_dp')

        phys_moisture.models.regenerate()

        # Single-phase conduit conductances of water liquid and vapour
        conduit_model = OpenPNM.Physics.models.multiphase.conduit_conductance
        tmodel = bpnm.Physics.models.diffusive_conductance.tdiff_conduit

        phys_water.models.add(propname='throat.conduit_hydraulic_conductance',
                              model=tmodel,
                              throat_conductance='throat.hydraulic_conductance',
                              pore_conductance='throat.hydraulic_conductance_pore')
        phys_water.models.add(model=conduit_model,
               propname='throat.conduit_hydraulic_conductance_wp',
               throat_conductance='throat.conduit_hydraulic_conductance',
               throat_occupancy='throat.occupancy_wp',
               pore_occupancy='pore.occupancy_wp', factor=1e-10)


        phys_vapour.models.add(propname='throat.conduit_diffusive_conductance',
                              model=tmodel,
                              throat_conductance='throat.diffusive_conductance_wp',
                              pore_conductance='throat.diffusive_conductance_wp_pore')
        phys_vapour.models.add(propname='throat.conduit_diffusive_conductance_wp',
                               model=conduit_model,
                               throat_conductance='throat.conduit_diffusive_conductance',
                               throat_occupancy='throat.occupancy_wp',
                               pore_occupancy='pore.occupancy_wp', factor=1e-10)

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp + dPc
#        pc2_wp = Pc_step_wp - dPc
        pc1_dp = Pc_step_dp + dPc
#        pc2_dp = Pc_step_dp - dPc

        for bound_increment in range(len(bounds)):
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            alg_flow_dp = pab.MoistureFlow(name='alg_flow_dp', network=pn,
                                           phase=moisture)
            alg_water_wp = pab.WaterFlow(name='alg_water_wp',
                                         network=pn, phase=water)
            alg_vapour_wp = pab.WaterFlow(name='alg_vapour_wp',
                                          network=pn, phase=vapour)

            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            Algs = [[alg_flow_wp, pc1_wp], [alg_flow_dp, pc1_dp],
                    [alg_water_wp, pc1_wp], [alg_vapour_wp, pc1_wp]]
            # BC1
            for alg, pc1 in Algs:
                alg.set_boundary_conditions(bctype='Dirichlet',
                                            bcvalue=pc1, pores=BC1_pores)

            # BC2
            for alg, pc1 in Algs:
                alg.set_boundary_conditions(bctype='Dirichlet',
                                            bcvalue=pc1-2*dPc, pores=BC2_pores)

            # run algorithms with proper conduit conductance
            alg_flow_wp.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')
            alg_flow_dp.run(conductance='conduit_conductance_dp',
                            quantity='pressure_dp')
            alg_water_wp.run(conductance='conduit_hydraulic_conductance_wp',
                            quantity='pressure_wp')
            alg_vapour_wp.run(conductance='conduit_diffusive_conductance_wp',
                            quantity='pressure_wp')

#            alg_flow_dp.return_results()
            # calc effective permeabilities [s]
            eff_permeability_moisture_wp = alg_flow_wp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            eff_permeability_moisture_dp = alg_flow_dp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_dp'])
            eff_permeability_water_wp = alg_water_wp.calc_eff_permeability(
                conductance=phys_water['throat.conduit_hydraulic_conductance_wp'])
            eff_permeability_vapour_wp = alg_vapour_wp.calc_eff_permeability(
                conductance=phys_vapour['throat.conduit_diffusive_conductance_wp'])

            # append permeability & flow values to the lists
            eff_perm_moisture_wp[str(bound_increment)].append(
                eff_permeability_moisture_wp)
            eff_perm_moisture_dp[str(bound_increment)].append(
                eff_permeability_moisture_dp)
            eff_perm_water_wp[str(bound_increment)].append(
                eff_permeability_water_wp)
            eff_perm_vapour_wp[str(bound_increment)].append(
                eff_permeability_vapour_wp)

            ctrl.purge_object(alg_flow_wp)
            ctrl.purge_object(alg_flow_dp)
            ctrl.purge_object(alg_water_wp)
            ctrl.purge_object(alg_vapour_wp)

#==============================================================================
#             #% Verification: compare water occupancy
#             # --------------------------------------
#             Pc_p = alg_flow_wp['pore.moisture_pressure_wp']     # Pc result
#             connected_pores = pn['throat.conns']
#             Pc_connected_pore = [[Pc_p[pair]] for pair in connected_pores]
#             Pc_connected_pore = sp.array(Pc_connected_pore).reshape(
#                                         (sp.shape(connected_pores)))
#             Pc_t_result = sp.amin(Pc_connected_pore, axis=1)
# #            Pc_t_result = pn.interpolate_data(data=Pc_p)  # mean of 2 pores
#
#             Tinvaded1 = water['throat.occupancy_wp']
#             Tinvaded2 = sp.float64(alg_wp._t_cap<=Pc_step_wp)
#             Tinvaded3 = sp.float64(alg_wp._t_cap<=Pc_t_result)
#             diff12 = sp.where(Tinvaded1!=Tinvaded2)
#             diff23 = sp.where(Tinvaded3!=Tinvaded2)
#             if sp.size(diff12):
#                 print('Different12 water distribution at: ',
#                       diff12)
#                 print('Pc throat: ', alg_wp._t_cap[diff12])
#                 print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
#                 print('Pc step throat: ', Pc_t_result[diff12])
#                 print('Pc step conn pores: ', Pc_connected_pore[diff12])
#
#             if sp.size(diff23):
#                 print('Different23 water distribution at: ',
#                       diff23)
#                 print('Pc throat: ', alg_wp._t_cap[diff23])
#                 print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
#                 print('Pc step throat: ', Pc_t_result[diff23])
#                 print('Pc step conn pores: ', Pc_connected_pore[diff23])
#             Ax = alg_flow_wp.A.dot(alg_flow_wp.X)
#             b = alg_flow_wp.b.reshape(sp.shape(Ax))
#             res = Ax - b
# #            norm_res = sp.linalg.norm(res, ord=2)
# #            print('Residual 2-norm: ', norm_res, '\n')
# #            max_norm_res = sp.amax(max_norm_res, norm_res)
#==============================================================================

    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=eff_perm_moisture_wp['0'])
    alg_flow_dp.store_result(Pc=Pc_dp, sat=sat_dp, sat_surf=sat_dp_surf,
                             sat_moist=sat_dp_moist,
                             w_sat=w_sat, k=eff_perm_moisture_dp['0'])
    alg_water_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=eff_perm_water_wp['0'],
                              span_water=span_water, span_air=span_air)
    alg_vapour_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=eff_perm_vapour_wp['0'],
                              span_water=span_water, span_air=span_air)

#    alg_flow_wp.calc_abs_permeability()
#    alg_flow_wp.calc_mD_permeability()
#    alg_flow_dp.calc_abs_permeability()
#    alg_flow_dp.calc_mD_permeability()
#    alg_water_wp.calc_abs_permeability()
#    alg_water_wp.calc_mD_permeability()

    if plot:
        bpnm.Postprocessing.Plots.hysteresis(alg_flow_wp, alg_flow_dp)
    return (alg_flow_wp, alg_flow_dp, alg_water_wp, alg_vapour_wp)


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
