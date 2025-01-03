# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 12:10:39 2015

@author: islah

This is a module which consists of routine functions for moisture storage and transfer estimations for topological network model

--> use throat.porelengths in conductivity calculation

"""
import bwfpnm as bpnm
import scipy as sp

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
                                  pores=ps, throats=ts,
                                  dynamic_data=True, name=physvapourname)
    phys_moisture = PhysMoistClass(network=NetInstance, phase=MoistInstance,
                                   pores=ps, throats=ts,
                                   dynamic_data=True, name=physmoistname)

    return (phys_water, phys_vapour, phys_moisture)


#%% Create algorithm: Wetting & Drying Percolation -> only use element's Pc_wp
def algorithm(NetIns, WaterIns, VapourIns, AdsorpClass=None, DesorpClass=None,
              npts=None, inv_points=None, air_entrapment=False,
              inv_sites=None, **kwargs):

    if AdsorpClass is None:
        AdsorpClass = bpnm.Algorithms.WettingPercolation
    alg_wp = AdsorpClass(network=NetIns,
                         invading_phase=WaterIns, defending_phase=VapourIns,
                         name='WettingPercolation', **kwargs)
    if inv_sites is None:
        inv_sites = NetIns['pore.inlet'] + NetIns['pore.outlet']    # bool
    alg_wp.run(inlets=inv_sites, npts=npts, inv_points=inv_points,
               access_limited=air_entrapment)

    return alg_wp


#%% create moisture content (w) array
def moist_content(geo, alg_wp, water_density,
                  v_mat=None, porosity=None, **kwargs):

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

    alg_wp['pore.inv_w'] = alg_wp['pore.inv_sat']*w_sat

    try:
        alg_wp['throat.inv_w'] = alg_wp['throat.inv_sat']*w_sat
    except:
        pass

    return (w_sat, porosity)


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
def permeability(pn, alg_wp, water, vapour, moisture,
                 phys_vapour, phys_moisture, w_sat,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc =1,
                 **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat_wp = []
    sat_wp_surf = []
    sat_wp_moist = []
    max_norm_res = 0
    eff_perm_moisture_wp = {'0': [], '1': [], '2': []}


    p_volumes = pn['pore.volume']
    t_volumes = pn['throat.volume']
    volume_total = sum(p_volumes) + sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
    Pc_wp = sp.around(Pc_wp, 3)

    for Pc_step_wp in Pc_wp:
        alg_wp.return_results(Pc=Pc_step_wp, occupancy='occupancy_wp')

        p_occ_wp = water['pore.occupancy_wp']
        t_occ_wp = water['throat.occupancy_wp']

        volume_p_wp = sum(p_occ_wp*p_volumes)
        volume_t_wp = sum(t_occ_wp*t_volumes)

        saturation_wp = (volume_p_wp + volume_t_wp)/(volume_total)
        sat_wp.append(saturation_wp)

        if printstatus:
            print('WP_saturation: %.3f' % saturation_wp,
                  '\t Pc: %.3f' %Pc_step_wp)
            print('WP_volume: ', sum(p_occ_wp),
                  '\t throat: ', sum(t_occ_wp))
            print('WP_water occupancy: ', p_occ_wp,
                  '\t throat: ', t_occ_wp)

        # Surface adsorption: p/t.surface_thickness & surface_volume
        if surface_ad:
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_wp)

            ## Inner circular tube
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume)
            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume)

            volume_p_wp = sum(phys_vapour['pore.surface_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])

            sat_surf_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)

            sat_wp[-1] += sat_surf_wp   # update total saturations

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

            volume_p_wp = sum(phys_vapour['pore.moist_volume_wp'])
            volume_t_wp = sum(phys_vapour['throat.moist_volume_wp'])

            sat_moist_wp = (volume_p_wp + volume_t_wp)/(volume_total)
            sat_wp_moist.append(sat_moist_wp)

            sat_wp[-1] += sat_moist_wp   # update total saturations

            print('moist vol: ', phys_vapour['throat.moist_volume_wp'])
            print('moist sat: ',
                  phys_vapour['throat.moist_volume_wp']/t_volumes)

        # Update vapour permeability for ALL pores & throats for wp & dp
        # g = f(Pv(rh(pc)))
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp_pore',
                                model=pm.diffusive_conductance.tbulk_diffusion_pore,
                                pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.regenerate()

        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance_pore,
                                 throat_occupancy='throat.occupancy_wp',
                                 pore_occupancy='pore.occupancy_wp',
                                 pdiffusive_conductance='throat.diffusive_conductance_wp_pore',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp')

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp + dPc
        pc2_wp = Pc_step_wp - dPc

        for bound_increment in range(len(bounds)):
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture,
                                           conductance='conduit_conductance_wp',
                                           quantity='pressure_wp')

            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # BC1
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_wp,
                                                pores=BC1_pores)

            # BC2
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_wp,
                                                pores=BC2_pores)

            # run algorithms with proper conduit conductance
            alg_flow_wp.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')

            # calc effective permeabilities [s]
            eff_permeability_moisture_wp = alg_flow_wp.calc_eff_permeability(
                conductance=phys_moisture['throat.conduit_conductance_wp'])

            # append permeability & flow values to the lists
            eff_perm_moisture_wp[str(bound_increment)].append(
                eff_permeability_moisture_wp)

            ctrl.purge_object(alg_flow_wp)

            #% Verification: compare water occupancy
            # --------------------------------------
            Pc_p = alg_flow_wp['pore.moisture_pressure_wp']     # Pc result
            connected_pores = pn['throat.conns']
            Pc_connected_pore = [[Pc_p[pair]] for pair in connected_pores]
            Pc_connected_pore = sp.array(Pc_connected_pore).reshape(
                                        (sp.shape(connected_pores)))
            Pc_t_result = sp.amin(Pc_connected_pore, axis=1)
#            Pc_t_result = pn.interpolate_data(data=Pc_p)  # mean of 2 pores

            Tinvaded1 = water['throat.occupancy_wp']
            Tinvaded2 = sp.float64(alg_wp._t_cap<=Pc_step_wp)
            Tinvaded3 = sp.float64(alg_wp._t_cap<=Pc_t_result)
            diff12 = sp.where(Tinvaded1!=Tinvaded2)
            diff23 = sp.where(Tinvaded3!=Tinvaded2)
            if sp.size(diff12):
                print('Different12 water distribution at: ',
                      diff12)
                print('Pc throat: ', alg_wp._t_cap[diff12])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: ', Pc_t_result[diff12])
                print('Pc step conn pores: ', Pc_connected_pore[diff12])

            if sp.size(diff23):
                print('Different23 water distribution at: ',
                      diff23)
                print('Pc throat: ', alg_wp._t_cap[diff23])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: ', Pc_t_result[diff23])
                print('Pc step conn pores: ', Pc_connected_pore[diff23])
            Ax = alg_flow_wp.A.dot(alg_flow_wp.X)
            b = alg_flow_wp.b.reshape(sp.shape(Ax))
            res = Ax - b
#            norm_res = sp.linalg.norm(res, ord=2)
#            print('Residual 2-norm: ', norm_res, '\n')
#            max_norm_res = sp.amax(max_norm_res, norm_res)

    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, k=eff_perm_moisture_wp['0'])

    alg_flow_wp.calc_abs_permeability()
    alg_flow_wp.calc_mD_permeability()

    if plot:
        bpnm.Postprocessing.Plots.plot_wp(alg_flow_wp)
    return alg_flow_wp


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
