# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 20:32:12 2015

@author: islah
"""
import bwfpnm as bpnm
import scipy as sp
from OpenPNM.Base import logging
from numpy.linalg import cond as np_cond
from multiprocessing import Pool
logger = logging.getLogger(__name__)

ctrl = bpnm.Base.Controller()
ctrl.loglevel = 40

def permeability_wp(pn, alg_wp, water, vapour, moisture,
                    phys_vapour, phys_moisture, w_sat,
                    num_seq=10, knudsen=True, plot=True, printstatus=False,
                    surface_ad=False, moist_volume=False, dPc = 1,
                    inv_points=None, n_cpu=4, **kwargs):

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

    args = [{'Pc_step_wp': Pc_step} for Pc_step in Pc_wp]

    pool = Pool(n_cpu)
    pool.map(permeability_main, args)
    pool.close()
    pool.join()

    for Pc_step_wp in Pc_wp:
        permeability_main()


    alg_flow_wp.store_result(Pc=Pc_wp, sat=sat_wp, sat_surf=sat_wp_surf,
                             sat_moist=sat_wp_moist,
                             w_sat=w_sat, g=eff_conduct_wp['0'])
    print('Total residual 2-norm: ', max_norm_res, '\n')

    return alg_flow_wp


def permeability_main(Pc_step_wp, t_volumes, volume_total, sat_wp, surface_ad,
                      pm, sat_wp_surf, moist_volume, sat_wp_moist, print_status,
                      knudsen, dPc, ):
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
        print('Pc step: ', Pc_step_wp)
        print('Relative error: ', rel_error)
        print('Residual 2-norm: ', norm_res, '\n')

#            print('Condition number: ', np_cond(alg_flow_wp.A.todense()))
        max_norm_res = sp.amax(max_norm_res, norm_res)

#    if saturation_wp==1:
#        # STOP calculation when SATuration = 1
#        n = num_seq - len(sat_wp)
#        sat_wp.extend([saturation_wp]*n)
#        eff_conduct_wp[str(bound_increment)].extend(
#            [eff_conduct_moisture_wp]*n)
#        if surface_ad:
#            sat_wp_surf.extend([sat_surf_wp]*n)
#        if moist_volume:
#            sat_wp_moist.extend([sat_moist_wp]*n)
#        break

# %% Main
if __name__ == '__main__':
    (pn, geo, water, vapour, moisture, phys_water, phys_vapour, phys_moisture,
     alg_wp, alg_dp, w_sat, porosity) = moisture_storage(plot=False)
    alg_flow_wp, alg_flow_dp = permeability(pn, alg_wp, alg_dp,
                                            water, vapour, moisture,
                                            phys_vapour, phys_moisture,
                                            num_seq=10, w_sat=w_sat)
