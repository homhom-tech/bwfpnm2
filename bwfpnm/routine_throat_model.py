# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:04:26 2015

@author: islah
"""
import bwfpnm as bpnm
import scipy as sp
from OpenPNM.Base import logging
from numpy.linalg import cond as np_cond
logger = logging.getLogger(__name__)

ctrl = bpnm.Base.Controller()
ctrl.loglevel = 40


# %% Network
def network(NetClass=None, netname='net', **kwargs):

    # Create network object
    if NetClass is None:
        NetClass = bpnm.Network.Cubic

    pn = NetClass(name=netname, **kwargs)
    return pn


# %% Geometry
def geometry(NetInstance, GeoClass=None, geoname='geo',
             pores='all', throats='all', **kwargs):

    ps = NetInstance.pores(pores)
    ts = NetInstance.throats(throats)

    if GeoClass is None:
        try:
            GeoClass = bpnm.Geometry.Throat
        except:
            GeoClass = bpnm.Geometry.Test_PoreThroat

    geo = GeoClass(network=NetInstance, pores=ps, throats=ts,
                   name=geoname, **kwargs)

    return geo


# %% Create phases: water + vapour = moisture
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


# %% Create physics
def physic(NetInstance, WaterInstance, VapourInstance, MoistInstance,
           pores='all', throats='all',
           PhysWaterClass=None, PhysVapourClass=None, PhysMoistClass=None,
           physwatername='physics_water', physvapourname='physics_vapour',
           physmoistname='physics_moisture', **kwargs):

    if PhysWaterClass is None:
        PhysWaterClass = bpnm.Physics.Standard
    if PhysVapourClass is None:
        PhysVapourClass = bpnm.Physics.GenericPhysics
    if PhysMoistClass is None:
        PhysMoistClass = bpnm.Physics.GenericPhysics

    ps = NetInstance.pores(pores)
    ts = NetInstance.throats(throats)
    phys_water = PhysWaterClass(network=NetInstance, phase=WaterInstance,
                                pores=ps, throats=ts,
                                dynamic_data=True, name=physwatername)
    phys_vapour = PhysVapourClass(network=NetInstance, phase=VapourInstance,
                                  pores=ps, throats=ts,
                                  dynamic_data=True, name=physvapourname)
    phys_moisture = PhysMoistClass(network=NetInstance, phase=MoistInstance,
                                   pores=ps, throats=ts,
                                   dynamic_data=True, name=physmoistname)

    return (phys_water, phys_vapour, phys_moisture)


# %% Create algorithm: Wetting & Drying Percolation -> only use element's Pc_wp
def algorithm(NetIns, WaterIns, VapourIns, AdsorpClass=None, DesorpClass=None,
              npts=None, inv_points=None, dp=True, **kwargs):

    if AdsorpClass is None:
        AdsorpClass = bpnm.Algorithms.WettingPercolation
    alg_wp = AdsorpClass(network=NetIns,
                         invading_phase=WaterIns, defending_phase=VapourIns,
                         name='WettingPercolation', **kwargs)
    inv_sites = NetIns['pore.inlet'] + NetIns['pore.outlet']    # bool
    alg_wp.run(inlets=None, npts=npts, inv_points=inv_points)

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


# %% create moisture content (w) array
def moist_content(geo, alg_wp, alg_dp, water_density,
                  v_mat=None, porosity=None, **kwargs):

    try:
        v_pore = sp.sum(geo['pore.volume']) + sp.sum(geo['throat.volume'])
    except:
        v_pore = sp.sum(geo['pore.volume'])

    if v_mat is not None:
        porosity = v_pore/v_mat

    try:
        w_sat = porosity*water_density
    except:
        print('error: either volume of bulk material of porosity is required!')
        return

    alg_wp['pore.inv_w'] = alg_wp['pore.inv_sat']*w_sat
    alg_dp['pore.inv_w'] = alg_dp['pore.inv_sat']*w_sat
    try:
        alg_wp['throat.inv_w'] = alg_wp['throat.inv_sat']*w_sat
        alg_dp['throat.inv_w'] = alg_dp['throat.inv_sat']*w_sat
    except:
        pass

    return (alg_wp, alg_dp, w_sat, porosity)


# %% Postprocessing: Plot
def plot_moist_retention(alg_wp, sat=None, Pc=None, **kwargs):

    if sat is not None:
        alg_wp.return_results(sat=sat)
    elif Pc is not None:
        alg_wp.return_results(Pc=Pc)
    else:
        alg_wp.return_results(sat=.5)
#     bpnm.Postprocessing.Plots.moisture_retention(alg_wp, w_sat)
    bpnm.Postprocessing.Plots.wetting_curves(alg_wp)


# %% Import for Paraview
def create_vtk(NetInstance, listPhases, filename='wetting.vtp', **kwargs):
    import bwfpnm.Utilities.IO as io
    io.VTK.save(network=NetInstance, filename=filename, phases=listPhases)


# %% Calculate permeability for each Pc value
def permeability(pn, alg_wp, alg_dp, water, vapour, moisture,
                 phys_vapour, phys_moisture, w_sat,
                 num_seq=10, knudsen=True, plot=True, printstatus=False,
                 surface_ad=False, moist_volume=False, dPc = 1,
                 **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat_wp, sat_dp = [], []
    sat_wp_surf, sat_dp_surf = [], []
    sat_wp_moist, sat_dp_moist = [], []

    eff_conduct_wp = {'0': [], '1': [], '2': []}
    eff_conduct_dp = {'0': [], '1': [], '2': []}

    t_volumes = pn['throat.volume']
    volume_total = sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
    Pc_wp = sp.around(Pc_wp, 3)
    Pc_dp = Pc_wp[::-1]

    for Pc_step_wp, Pc_step_dp in list(zip(Pc_wp, Pc_dp)):
        alg_wp.return_results(Pc=Pc_step_wp, occupancy='occupancy_wp')
        alg_dp.return_results(Pc=Pc_step_dp, occupancy='occupancy_dp')

        t_occ_wp = water['throat.occupancy_wp']
        t_occ_dp = water['throat.occupancy_dp']

        volume_t_wp = sum(t_occ_wp*t_volumes)
        volume_t_dp = sum(t_occ_dp*t_volumes)

        saturation_wp = (volume_t_wp)/(volume_total)
        saturation_dp = (volume_t_dp)/(volume_total)
        sat_wp.append(saturation_wp)
        sat_dp.append(saturation_dp)

        # Surface adsorption: p/t.surface_thickness & surface_volume
        if surface_ad:  # only applied to dry pores
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.surface_thickness_dp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_dp')
            phys_vapour.models.add(propname='throat.surface_thickness_dp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_dp')

            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume,
                                   film_thickness='throat.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume,
                                   film_thickness='pore.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_dp',
                                   model=pm.surface_adsorption.pvolume,
                                   pore_occupancy='pore.occupancy_dp',
                                   film_thickness='pore.surface_thickness_dp')
            phys_vapour.models.add(propname='throat.surface_volume_dp',
                                   model=pm.surface_adsorption.tvolume,
                                   throat_occupancy='throat.occupancy_dp',
                                   film_thickness='throat.surface_thickness_dp')

            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])
            volume_t_dp = sum(phys_vapour['throat.surface_volume_dp'])

            sat_surf_wp = (volume_t_wp)/(volume_total)
            sat_surf_dp = (volume_t_dp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)
            sat_dp_surf.append(sat_surf_dp)

            sat_wp[-1] += sat_surf_wp   # update total saturations
            sat_dp[-1] += sat_surf_dp

            print('tthickness wp: ',
                  phys_vapour['throat.surface_thickness_wp'])
            print('tsat wp: ',
                  phys_vapour['throat.surface_volume_wp']/t_volumes)

        if moist_volume:
            phys_vapour.models.add(propname='pore.moist_volume_wp',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='throat.moist_volume_wp',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.moist_volume_dp',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step_dp,
                                   pore_occupancy='pore.occupancy_dp')
            phys_vapour.models.add(propname='throat.moist_volume_dp',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step_dp,
                                   throat_occupancy='throat.occupancy_dp')

            volume_t_wp = sum(phys_vapour['throat.moist_volume_wp'])
            volume_t_dp = sum(phys_vapour['throat.moist_volume_dp'])

            sat_moist_wp = (volume_t_wp)/(volume_total)
            sat_moist_dp = (volume_t_dp)/(volume_total)
            sat_wp_moist.append(sat_moist_wp)
            sat_dp_moist.append(sat_moist_dp)

            sat_wp[-1] += sat_moist_wp   # update total saturations
            sat_dp[-1] += sat_moist_dp

            print('moist vol: ', phys_vapour['throat.moist_volume_wp'])
            print('moist sat: ',
                  phys_vapour['throat.moist_volume_wp']/t_volumes)

        if printstatus:
#            saturation = saturation_wp + sat_surf_wp + sat_moist_wp
#            print('WP total saturation: %.3f' % saturation)
            print('WP_saturation: %.3f' % saturation_wp,
                  '\t Pc: %.3f' % sp.log10(-Pc_step_wp))
            print('DP_saturation: %.3f' % saturation_dp,
                  '\t Pc: %.3f' % sp.log10(-Pc_step_dp))
            print('WP_volume_throat: ', sum(t_occ_wp))
            print('DP_volume_ throat: ', sum(t_occ_dp))

        # Update vapour permeability for all pores & throats for wp & dp
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.add(propname='pore.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.pbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.add(propname='throat.diffusive_conductance_dp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_dp, knudsen=knudsen)
        phys_vapour.models.add(propname='pore.diffusive_conductance_dp',
                               model=pm.diffusive_conductance.pbulk_diffusion,
                               pc=Pc_step_dp, knudsen=knudsen)
        phys_vapour.models.regenerate()
        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance,
                                 throat_occupancy='throat.occupancy_wp',
                                 pore_occupancy='pore.occupancy_wp',
                                 pdiffusive_conductance='pore.diffusive_conductance_wp',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp')
        phys_moisture.models.add(propname='throat.conduit_conductance_dp',
                                 model=pm.multiphase.mixed_conductance,
                                 throat_occupancy='throat.occupancy_dp',
                                 pore_occupancy='pore.occupancy_dp',
                                 pdiffusive_conductance='pore.diffusive_conductance_dp',
                                 tdiffusive_conductance='throat.diffusive_conductance_dp')

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp - dPc/2
        pc2_wp = Pc_step_wp + dPc/2
        pc1_dp = Pc_step_dp - dPc/2
        pc2_dp = Pc_step_dp + dPc/2

        for bound_increment in range(len(bounds)):
            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # Flow for Wetting Percolation
            # ----------------------------
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            # BC1
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_wp,
                                                pores=BC1_pores)
            # BC2
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_wp,
                                                pores=BC2_pores)
            # run algorithms (pressure distribution, flow rate) with proper conduit conductance
            alg_flow_wp.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')
            # calc effective conductance: G/dP
            eff_conduct_moisture_wp = alg_flow_wp.calc_eff_conduct_conceptual(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            # append permeability & flow values to the lists
            eff_conduct_wp[str(bound_increment)].append(
                eff_conduct_moisture_wp)
            ctrl.purge_object(alg_flow_wp)

            # Flow for Drying Percolation
            # ---------------------------
            alg_flow_dp = pab.MoistureFlow(name='alg_flow_dp', network=pn,
                                           phase=moisture)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_dp,
                                                pores=BC1_pores)
            alg_flow_dp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_dp,
                                                pores=BC2_pores)
            alg_flow_dp.run(conductance='conduit_conductance_dp',
                            quantity='pressure_dp')
            eff_conduct_moisture_dp = alg_flow_dp.calc_eff_conduct_conceptual(
                conductance=phys_moisture['throat.conduit_conductance_dp'])
            eff_conduct_dp[str(bound_increment)].append(
                eff_conduct_moisture_dp)
            ctrl.purge_object(alg_flow_dp)

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
                print('Different12 water distribution at: \n',
                      diff12)
                print('Pc throat: ', alg_wp._t_cap[diff12])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: \n', Pc_t_result[diff12])
                print('Pc step conn pores: \n', Pc_connected_pore[diff12], '\n')

            if sp.size(diff23):
                print('Different23 water distribution at: \n',
                      diff23)
                print('Pc throat: \n', alg_wp._t_cap[diff23])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: \n', Pc_t_result[diff23])
                print('Pc step conn pores: \n', Pc_connected_pore[diff23], '\n')
        if saturation_wp==1:
            # STOP calculation when SATuration = 1
            n = num_seq - len(sat_wp)
            sat_wp.extend([saturation_wp]*n)
            eff_conduct_wp[str(bound_increment)].extend(
                [eff_conduct_moisture_wp]*n)
            if surface_ad:
                sat_wp_surf.extend([sat_surf_wp]*n)
            if moist_volume:
                sat_wp_moist.extend([sat_moist_wp]*n)
            break

    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, g=eff_conduct_wp['0'])
    alg_flow_dp.store_result(Pc=Pc_dp, sat=sat_dp, sat_surf=sat_dp_surf,
                             sat_moist=sat_dp_moist,
                             w_sat=w_sat, g=eff_conduct_dp['0'])

    if plot:
        bpnm.Postprocessing.Plots.hysteresis(alg_flow_wp, alg_flow_dp)
    return (alg_flow_wp, alg_flow_dp)


def permeability_wp(pn, alg_wp, water, vapour, moisture,
                    phys_vapour, phys_moisture, w_sat,
                    num_seq=10, knudsen=True, plot=True, printstatus=False,
                    surface_ad=False, moist_volume=False, dPc = 1,
                    inv_points=None, **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat_wp = []
    sat_wp_surf = []
    sat_wp_moist = []
    eff_conduct_wp = {'0': [], '1': [], '2': []}
    max_norm_res = 0

    t_volumes = pn['throat.volume']
    volume_total = sum(t_volumes)

    Pc_wp = inv_points
    if inv_points is None:
        lr = sp.arange(-10, -1)
        r = sp.power(10, lr)
        pc = -2*water['pore.surface_tension'][0]/r

        Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
        Pc_wp = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
        Pc_wp = sp.around(Pc_wp, 3)
    else:
        num_seq = sp.size(inv_points)

    for Pc_step_wp in Pc_wp:
        alg_wp.return_results(Pc=Pc_step_wp, occupancy='occupancy_wp')

        t_occ_wp = water['throat.occupancy_wp']
        volume_t_wp = sum(t_occ_wp*t_volumes)

        saturation_wp = (volume_t_wp)/(volume_total)
        sat_wp.append(saturation_wp)

        # Surface adsorption: p/t.surface_thickness & surface_volume
        if surface_ad:  # only applied to dry pores
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step_wp)

            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume,
                                   film_thickness='throat.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume,
                                   film_thickness='pore.surface_thickness_wp')

            volume_t_wp = sum(phys_vapour['throat.surface_volume_wp'])

            sat_surf_wp = (volume_t_wp)/(volume_total)
            sat_wp_surf.append(sat_surf_wp)

            print('tthickness wp: ',
                  phys_vapour['throat.surface_thickness_wp'])
            print('tsat wp: ',
                  phys_vapour['throat.surface_volume_wp']/t_volumes)

        if moist_volume:
            phys_vapour.models.add(propname='throat.moist_volume_wp',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step_wp)
            phys_vapour.models.add(propname='pore.moist_volume_wp',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step_wp)

            volume_t_wp = sum(phys_vapour['throat.moist_volume_wp'])

            sat_moist_wp = (volume_t_wp)/(volume_total)
            sat_wp_moist.append(sat_moist_wp)

            print('moist vol: ', phys_vapour['throat.moist_volume_wp'])
            print('moist sat: ',
                  phys_vapour['throat.moist_volume_wp']/t_volumes)

        if printstatus:
#            saturation = saturation_wp + sat_surf_wp + sat_moist_wp
#            print('WP total saturation: %.3f' % saturation)
            print('WP_saturation: %.3f' % saturation_wp,
                  '\t Pc: %.3f' % sp.log10(-Pc_step_wp))
            print('WP_volume_throat: ', sum(t_occ_wp))

        # Update vapour permeability for all pores & throats for wp & dp
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)
        phys_vapour.models.add(propname='pore.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.pbulk_diffusion,
                               pc=Pc_step_wp, knudsen=knudsen)

        phys_vapour.models.regenerate()
        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance,
                                 throat_occupancy='throat.occupancy_wp',
                                 pore_occupancy='pore.occupancy_wp',
                                 pdiffusive_conductance='pore.diffusive_conductance_wp',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp')

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1_wp = Pc_step_wp - dPc/2
        pc2_wp = Pc_step_wp + dPc/2

        for bound_increment in range(len(bounds)):
            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # Flow for Wetting Percolation
            # ----------------------------
            alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            # BC1
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1_wp,
                                                pores=BC1_pores)
            # BC2
            alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2_wp,
                                                pores=BC2_pores)
            # run algorithms (pressure distribution, flow rate) with proper conduit conductance
            alg_flow_wp.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')
            # calc effective conductance: G/dP
            eff_conduct_moisture_wp = alg_flow_wp.calc_eff_conduct_conceptual(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            # append permeability & flow values to the lists
            eff_conduct_wp[str(bound_increment)].append(
                eff_conduct_moisture_wp)
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
                # FIXED by supplying the same 'inv_points' argument to percolation and permeability algorithms
                print('Different12 water distribution at: \n',
                      diff12)
                print('Pc throat: ', alg_wp._t_cap[diff12])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: \n', Pc_t_result[diff12])
                print('Pc step conn pores: \n', Pc_connected_pore[diff12], '\n')

            if sp.size(diff23):
                print('Different23 water distribution at: \n',
                      diff23)
                print('Pc throat: \n', alg_wp._t_cap[diff23])
                print('Pc step wp: ', Pc_step_wp, 'lPc: ', sp.log10(-Pc_step_wp))
                print('Pc step throat: \n', Pc_t_result[diff23])
                print('Pc step conn pores: \n', Pc_connected_pore[diff23], '\n')

            Ax = alg_flow_wp.A.dot(alg_flow_wp.X)
            b = alg_flow_wp.b.reshape(sp.shape(Ax))
            res = Ax - b
            norm_res = sp.linalg.norm(res, ord=2)
            rel_error = sp.linalg.norm(Ax)/sp.linalg.norm(b)
            print('Relative error: ', rel_error)
            print('Residual 2-norm: ', norm_res, '\n')

#            print('Condition number: ', np_cond(alg_flow_wp.A.todense()))
            max_norm_res = sp.amax(max_norm_res, norm_res)

        if saturation_wp==1:
            # STOP calculation when SATuration = 1
            n = num_seq - len(sat_wp)
            sat_wp.extend([saturation_wp]*n)
            eff_conduct_wp[str(bound_increment)].extend(
                [eff_conduct_moisture_wp]*n)
            if surface_ad:
                sat_wp_surf.extend([sat_surf_wp]*n)
            if moist_volume:
                sat_wp_moist.extend([sat_moist_wp]*n)
            break

    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, g=eff_conduct_wp['0'])
    print('Total residual 2-norm: ', max_norm_res, '\n')

    return alg_flow_wp


def permeability_dp(pn, alg, water, vapour, moisture,
                    phys_vapour, phys_moisture, w_sat,
                    num_seq=10, knudsen=True, plot=True, printstatus=False,
                    surface_ad=False, moist_volume=False, dPc = 1,
                    **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    sat = []
    sat_surf = []
    sat_moist = []
    eff_conduct = {'0': [], '1': [], '2': []}

    t_volumes = pn['throat.volume']
    volume_total = sum(t_volumes)

    lr = sp.arange(-10, -1)
    r = sp.power(10, lr)
    pc = -2*water['pore.surface_tension'][0]/r

    Pc_min, Pc_max = sp.amin(pc), sp.amax(pc)
    Pc = -sp.logspace(sp.log10(-Pc_min), sp.log10(-Pc_max), num_seq)
    Pc = sp.around(Pc, 3)

    for Pc_step in Pc:
        alg.return_results(Pc=Pc_step, occupancy='occupancy')

        t_occ = water['throat.occupancy_wp']
        volume_t = sum(t_occ*t_volumes)

        saturation = (volume_t)/(volume_total)
        sat.append(saturation)

        # Surface adsorption: p/t.surface_thickness & surface_volume
        if surface_ad:  # only applied to dry pores
            phys_vapour.models.add(propname='throat.surface_thickness_wp',
                                   model=pm.surface_adsorption.tstat_thickness,
                                   pc=Pc_step)
            phys_vapour.models.add(propname='pore.surface_thickness_wp',
                                   model=pm.surface_adsorption.pstat_thickness,
                                   pc=Pc_step)

            phys_vapour.models.add(propname='throat.surface_volume_wp',
                                   model=pm.surface_adsorption.tvolume,
                                   film_thickness='throat.surface_thickness_wp')
            phys_vapour.models.add(propname='pore.surface_volume_wp',
                                   model=pm.surface_adsorption.pvolume,
                                   film_thickness='pore.surface_thickness_wp')

            volume_t = sum(phys_vapour['throat.surface_volume_wp'])

            sat_surf = (volume_t)/(volume_total)
            sat_surf.append(sat_surf)

            print('tthickness wp: ',
                  phys_vapour['throat.surface_thickness_wp'])
            print('tsat wp: ',
                  phys_vapour['throat.surface_volume_wp']/t_volumes)

        if moist_volume:
            phys_vapour.models.add(propname='throat.moist_volume_wp',
                                   model=pm.volume_moisture.tvolume,
                                   pc=Pc_step)
            phys_vapour.models.add(propname='pore.moist_volume_wp',
                                   model=pm.volume_moisture.pvolume,
                                   pc=Pc_step)

            volume_t = sum(phys_vapour['throat.moist_volume_wp'])

            sat_moist = (volume_t)/(volume_total)
            sat_moist.append(sat_moist)

            print('moist vol: ', phys_vapour['throat.moist_volume_wp'])
            print('moist sat: ',
                  phys_vapour['throat.moist_volume_wp']/t_volumes)

        if printstatus:
#            saturation = saturation + sat_surf + sat_moist
#            print('WP total saturation: %.3f' % saturation)
            print('WP_saturation: %.3f' % saturation,
                  '\t Pc: %.3f' % sp.log10(-Pc_step))
            print('WP_volume_throat: ', sum(t_occ))

        # Update vapour permeability for all pores & throats for wp & dp
        phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.tbulk_diffusion,
                               pc=Pc_step, knudsen=knudsen)
        phys_vapour.models.add(propname='pore.diffusive_conductance_wp',
                               model=pm.diffusive_conductance.pbulk_diffusion,
                               pc=Pc_step, knudsen=knudsen)

        phys_vapour.models.regenerate()
        # Calculate conduit conductances as a function of water distribution
        phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                                 model=pm.multiphase.mixed_conductance,
                                 throat_occupancy='throat.occupancy_wp',
                                 pore_occupancy='pore.occupancy_wp',
                                 pdiffusive_conductance='pore.diffusive_conductance_wp',
                                 tdiffusive_conductance='throat.diffusive_conductance_wp')

        phys_moisture.models.regenerate()

        bounds = [['inlet', 'outlet']]
        pc1 = Pc_step - dPc/2
        pc2 = Pc_step + dPc/2

        for bound_increment in range(len(bounds)):
            BC1_pores = pn.pores(labels=bounds[bound_increment][0])
            BC2_pores = pn.pores(labels=bounds[bound_increment][1])

            # Flow for Wetting Percolation
            # ----------------------------
            alg_flow = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                           phase=moisture)
            # BC1
            alg_flow.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc1,
                                                pores=BC1_pores)
            # BC2
            alg_flow.set_boundary_conditions(bctype='Dirichlet',
                                                bcvalue=pc2,
                                                pores=BC2_pores)
            # run algorithms (pressure distribution, flow rate) with proper conduit conductance
            alg_flow.run(conductance='conduit_conductance_wp',
                            quantity='pressure_wp')
            # calc effective conductance: G/dP
            eff_conduct_moisture = alg_flow.calc_eff_conduct_conceptual(
                conductance=phys_moisture['throat.conduit_conductance_wp'])
            # append permeability & flow values to the lists
            eff_conduct[str(bound_increment)].append(
                eff_conduct_moisture)
            ctrl.purge_object(alg_flow)

            #% Verification: compare water occupancy
            # --------------------------------------
            Pc_p = alg_flow['pore.moisture_pressure_wp']     # Pc result
            connected_pores = pn['throat.conns']
            Pc_connected_pore = [[Pc_p[pair]] for pair in connected_pores]
            Pc_connected_pore = sp.array(Pc_connected_pore).reshape(
                                        (sp.shape(connected_pores)))
            Pc_t_result = sp.amin(Pc_connected_pore, axis=1)
#            Pc_t_result = pn.interpolate_data(data=Pc_p)  # mean of 2 pores

            Tinvaded1 = water['throat.occupancy_wp']
            Tinvaded2 = sp.float64(alg._t_cap<=Pc_step)
            Tinvaded3 = sp.float64(alg._t_cap<=Pc_t_result)
            diff12 = sp.where(Tinvaded1!=Tinvaded2)
            diff23 = sp.where(Tinvaded3!=Tinvaded2)
            if sp.size(diff12):
                print('Different12 water distribution at: \n',
                      diff12)
                print('Pc throat: ', alg._t_cap[diff12])
                print('Pc step wp: ', Pc_step, 'lPc: ', sp.log10(-Pc_step))
                print('Pc step throat: \n', Pc_t_result[diff12])
                print('Pc step conn pores: \n', Pc_connected_pore[diff12], '\n')

            if sp.size(diff23):
                print('Different23 water distribution at: \n',
                      diff23)
                print('Pc throat: \n', alg._t_cap[diff23])
                print('Pc step wp: ', Pc_step, 'lPc: ', sp.log10(-Pc_step))
                print('Pc step throat: \n', Pc_t_result[diff23])
                print('Pc step conn pores: \n', Pc_connected_pore[diff23], '\n')
        if saturation==1:
            # STOP calculation when SATuration = 1
            n = num_seq - len(sat)
            sat.extend([saturation]*n)
            eff_conduct[str(bound_increment)].extend(
                [eff_conduct_moisture]*n)
            if surface_ad:
                sat_surf.extend([sat_surf]*n)
            if moist_volume:
                sat_moist.extend([sat_moist]*n)
            break

    alg_flow.store_result(Pc=Pc, sat=sat, sat_surf=sat_surf,
                             sat_moist=sat_moist,
                             w_sat=w_sat, g=eff_conduct['0'])

    return alg_flow


def permeability_bwf(pn, alg_wp, alg_dp, alg_flow_wp, alg_flow_dp,
                     water, vapour, moisture,
                     phys_water, phys_vapour, phys_moisture, w_sat,
                     num_seq=10, knudsen=True, **kwargs):

    pm = bpnm.Physics.models
    pab = bpnm.Algorithms

    Pc_max = -10  # sp.amax(alg_wp._inv_points)
    alg_wp.return_results(Pc=Pc_max, occupancy='occupancy_wp')

    p_volumes = pn['pore.volume']
    t_volumes = pn['throat.volume']
    volume_total = sum(p_volumes) + sum(t_volumes)

    p_occ_wp = water['pore.occupancy_wp']
    t_occ_wp = water['throat.occupancy_wp']
    volume_p_wp = sum(p_occ_wp*p_volumes)
    volume_t_wp = sum(t_occ_wp*t_volumes)
    saturation_wp = (volume_p_wp + volume_t_wp)/(volume_total)

    phys_water.models.add(propname='throat.hydraulic_conductance_bwf',
                          model=pm.hydraulic_conductance.thagen_poiseuille)
    phys_water.models.add(propname='pore.hydraulic_conductance_bwf',
                          model=pm.hydraulic_conductance.phagen_poiseuille)
    phys_water.models.regenerate()
    phys_vapour.models.add(propname='throat.diffusive_conductance_wp',
                           model=pm.diffusive_conductance.tbulk_diffusion,
                           pc=Pc_max, knudsen=knudsen)
    phys_vapour.models.add(propname='pore.diffusive_conductance_wp',
                           model=pm.diffusive_conductance.pbulk_diffusion,
                           pc=Pc_max, knudsen=knudsen)
    phys_vapour.models.regenerate()
    phys_moisture.models.add(propname='throat.conduit_conductance_wp',
                             model=pm.multiphase.mixed_conductance,
                             throat_occupancy='throat.occupancy_wp',
                             pore_occupancy='pore.occupancy_wp',
                             pdiffusive_conductance='pore.diffusive_conductance_wp',
                             tdiffusive_conductance='throat.diffusive_conductance_wp',
                             phydraulic_conductance='pore.hydraulic_conductance_bwf',
                             thydraulic_conductance='throat.hydraulic_conductance_bwf')


    phys_moisture.models.regenerate()

    bounds = [['inlet', 'outlet']]
    pc1 = Pc_max - 1
    pc2 = Pc_max + 1

    for bound_increment in range(len(bounds)):
        alg_flow_wp = pab.MoistureFlow(name='alg_flow_wp', network=pn,
                                       phase=moisture)

        BC1_pores = pn.pores(labels=bounds[bound_increment][0])
        BC2_pores = pn.pores(labels=bounds[bound_increment][1])

        # BC1
        alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                            bcvalue=pc1,
                                            pores=BC1_pores)
        # BC2
        alg_flow_wp.set_boundary_conditions(bctype='Dirichlet',
                                            bcvalue=pc2,
                                            pores=BC2_pores)
        alg_flow_wp.run(conductance='conduit_conductance_wp',
                        quantity='pressure_wp')
        # calc effective permeabilities [s]
        alg_flow_wp.calc_eff_permeability(
            conductance=phys_moisture['throat.conduit_conductance_wp'])
        alg_flow_wp.calc_abs_permeability()
        alg_flow_wp.calc_mD_permeability()

        # calc effective flow
        alg_flow_wp.calc_eff_conduct_conceptual(
            conductance=phys_moisture['throat.conduit_conductance_wp'])
        ctrl.purge_object(alg_flow_wp)

    return (saturation_wp, alg_flow_wp._abs_mD_permeability)


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
