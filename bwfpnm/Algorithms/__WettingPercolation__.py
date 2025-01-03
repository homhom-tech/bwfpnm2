"""
===============================================================================
module __OrdinaryPercolation__: Ordinary Percolation Algorithm
===============================================================================

"""

import scipy as sp
#import numpy as np
import matplotlib.pyplot as plt
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)

class WettingPercolation(GenericAlgorithm):
    r"""
    Simulates a capillary drainage experiment by looping through a list of
    capillary pressures.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    invading_phase : OpenPNM Phase Object
        The phase to be forced into the network at increasingly high pressures

    defending_phase : OpenPNM Phase Object, optional
        The phase originally residing in the network prior to invasion.  This
        is only necessary so that the pressure at which the phase is drained
        can be stored on the phase.

    name : string, optional
        The name to assign to the Algorithm Object

    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> geo = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
    >>> phase1 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phase2 = OpenPNM.Phases.TestPhase(network=pn)
    >>> phys1 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase1,pores=pn.pores(),throats=pn.throats())
    >>> phys2 = OpenPNM.Physics.TestPhysics(network=pn, phase=phase2,pores=pn.pores(),throats=pn.throats())
    >>> OP = OpenPNM.Algorithms.OrdinaryPercolation(network=pn, name='OP',invading_phase=phase1, defending_phase=phase2)
    >>> OP.run(inlets=pn.pores('top'))
    >>> med_Pc = sp.median(OP['pore.inv_Pc'])
    >>> OP.update_results(med_Pc)
    >>> print(len(phase1.pores('occupancy'))) #should return '71' filled pores if everything is working normally
    71

    To run this algorithm, use 'setup()' to provide the necessary simulation
    """

    def __init__(self,invading_phase=None,defending_phase=None,**kwargs):
        r"""

        """
        super(WettingPercolation,self).__init__(**kwargs)
        self._phase_inv = invading_phase
        self._phase_def = defending_phase
        logger.debug("Create Wetting Percolation Algorithm Object")

    def run(self,
            inlets,
            npts=None,
            inv_points=None,
            capillary_pressure='capillary_pressure',
            access_limited=False,
            trapping=False,
            **kwargs):
        r'''
        Parameters
        ----------
        inlets : array_like
            The list of pores which are the injection sources

        npts : int, optional
            The number of pressure points to apply.  The list of pressures
            is logarithmically spaced between the lowest and highest throat
            entry pressures in the network.

        inv_points : array_like, optional
            A list of specific pressure points to apply.

        trapping : boolean
            Wetting phase that is cut-off from the outlets becomes immobile.
            If outlet pores have not been provided then this argument is
            ignored.

        Notes
        -----
        The 'inlet' pores are initially filled with invading fluid to start the
        simulation.  To avoid the capillary pressure curve showing a non-zero
        starting saturation at low pressures, it is necessary to apply boundary
        pores that have 0 volume, and set these as the inlets.


        '''
        logger.info("Run module of WettingPercolation is executed with inlets: {}, air entrapment: {}".format(bool(sp.any(inlets)), access_limited))
        # Parse params
        self._inv_sites = inlets
        self._npts = npts
        self._p_cap = capillary_pressure  # Name of throat entry pressure prop
        self._AL = access_limited
        self._TR = trapping

        #Create pore and throat conditions lists to store inv_val at which each is invaded
        self._p_inv = sp.zeros((self._net.num_pores(),), dtype=float)
        self._p_inv.fill(sp.inf)
        self._p_seq = sp.zeros_like(self._p_inv, dtype=int)
        self._t_inv = sp.zeros((self._net.num_throats(),), dtype=float)
        self._t_inv.fill(sp.inf)
        self._t_seq = sp.zeros_like(self._t_inv, dtype=int)

        #Determine the invasion pressures to apply
        try:
            self._t_cap = self._phase_inv['throat.'+self._p_cap]
            _t_p_cap = sp.r_[self._t_cap]
        except:
            logger.error('Capillary pressure not assigned to invading phase\
                '+self._phase_inv.name +', check for capillary pressure in\
                defending phase '+self._phase_def.name +' instead')
        try:
            self._tp_cap = self._phase_inv['pore.'+self._p_cap]
            _t_p_cap = sp.r_[_t_p_cap, self._tp_cap]
        except:
            pass
            logger.error('Capillary pressure not assigned to pore invading\
                phase ' + self._phase_def.name + ' nor to invading phase '
                + self._phase_inv.name)

#        _t_p_cap = sp.r_[self._t_cap, self._tp_cap]
        if self._npts is None:
            self._npts = sp.size(sp.unique(_t_p_cap))

        if inv_points is None:
            min_p = sp.amin(_t_p_cap)#*1.02  # nudge min_p down slightly
            max_p = sp.amax(_t_p_cap)#*0.98  # bump max_p up slightly
            if min_p == 0:
                min_p = sp.linspace(min_p, max_p, self._npts)[1]
            self._inv_points = -sp.logspace(sp.log10(-min_p),
                                            sp.log10(-max_p), self._npts)
            self._inv_points = sp.around(self._inv_points, 3)
        else:
            inv_points = sp.array(inv_points)
            inv_points.sort()
            inv_points[0] = inv_points[0]#*1.02
            inv_points[-1] = inv_points[-1]#*0.98
            self._inv_points = sp.around(inv_points, 3)
            self._npts = sp.size(inv_points)
        self._do_outer_iteration_stage()

    def _do_outer_iteration_stage(self):
        #Generate curve from points
        for inv_val in self._inv_points:
            #Apply one applied pressure and determine invaded pores
#            self._logger.info('Applying capillary pressure: '+str(inv_val))
            self._do_one_inner_iteration(inv_val)
        #Store results using networks' get/set method
        self['pore.inv_Pc'] = self._p_inv
        self['throat.inv_Pc'] = self._t_inv
        self._pt_inv = sp.around(sp.r_[self._t_inv, self._p_inv], 3)

        #Find invasion sequence values (to correspond with IP algorithm)
        self._p_seq = sp.searchsorted(sp.unique(self._p_inv),self._p_inv)
        self._t_seq = sp.searchsorted(sp.unique(self._t_inv),self._t_inv)
        self._p_seq_all = sp.searchsorted(sp.unique(self._pt_inv), self._p_inv)
        self._t_seq_all = sp.searchsorted(sp.unique(self._pt_inv), self._t_inv)
        # Pores with trapped air is numbered as 1e+9
        pseq_inf = sp.where(self._p_inv==sp.inf)
        tseq_inf = sp.where(self._t_inv==sp.inf)
        self._p_seq[pseq_inf] = 1e+9
        self._t_seq[tseq_inf] = 1e+9
        self._p_seq_all[pseq_inf] = 1e+9
        self._t_seq_all[tseq_inf] = 1e+9
        self['pore.inv_seq_all'] = self._p_seq_all
        self['throat.inv_seq_all'] = self._t_seq_all
        self['pore.inv_seq'] = self._p_seq
        self['throat.inv_seq'] = self._t_seq
        #Calculate Saturations
        pvol = self._net['pore.volume']
        tvol = self._net['throat.volume']
        v_total = sp.sum(pvol) + sp.sum(tvol)
        sat = 0.
        self['pore.inv_sat'] = 0.
        self['throat.inv_sat'] = 0.
        for i in range(self._npts):
            inv_pores = sp.where(self._p_seq_all==i)[0]
            inv_throats = sp.where(self._t_seq_all==i)[0]
            v_liquid = (sum(pvol[inv_pores]) + sum(tvol[inv_throats]))
            new_sat = v_liquid/v_total
            sat += new_sat
            self['pore.inv_sat'][inv_pores] = sat
            self['throat.inv_sat'][inv_throats] = sat

    def _do_one_inner_iteration(self, inv_val):
        r"""
        Determine which throats & pores are invaded at a given applied capillary pressure
        """
        Pinvaded = self._tp_cap <= inv_val
        Tinvaded = self._t_cap <= inv_val
        if self._AL:
            pmask0 = Pinvaded.copy()     # possible-to-wet pores

            #Finding all pores that CAN BE invaded at specified pressure inv_val
            clusters = self._net.find_clusters(~Tinvaded)

            #Find all pores with at least 1 dry throat
            Ts = self._net.throats()
            P12 = self._net.find_connected_pores(Ts)
            temp = P12[~Tinvaded]
            temp = sp.hstack((temp[:,0],temp[:,1]))
            Pinvaded[temp] = True

            #Add injection sites to Pinvaded
            Pinvaded[self._inv_sites] = True
            #Clean up clusters (not invaded = -1, invaded >=0)
            clusters = clusters*(Pinvaded) - (~Pinvaded)
            #Identify which clusters are connected to invasion sites
            inv_clusters = sp.unique(clusters[self._inv_sites])

            pmask1 = sp.in1d(clusters, inv_clusters)  # accessable pores
            pmask = pmask1*pmask0

            #Determine Pc_invaded for throats as well
            temp = self._net['throat.conns']
            # tmask: a throat is accessible if one of the pores is accessible
            tmask = (pmask1[temp[:, 0]] + pmask1[temp[:, 1]])*(Tinvaded)

        else:
            pmask = Pinvaded
            tmask = Tinvaded
        #Store result of invasion step
        self._p_inv[(self._p_inv==sp.inf)*(pmask)] = inv_val
        self._t_inv[(self._t_inv==sp.inf)*(tmask)] = inv_val

    def evaluate_trapping(self, p_outlets):
        r"""
        Finds trapped pores and throats after a full ordinary
        percolation drainage has been run

        Parameters
        ----------
        outlets : array_like
            A list of pores that define the wetting phase outlets.
            Disconnection from these outlets results in trapping.

        """
        self['pore.trapped'] = sp.zeros([self.Np, ], dtype=float)
        self['throat.trapped'] = sp.zeros([self.Nt, ], dtype=float)
        try:
            # Get points used in OP
            inv_points = sp.unique(self['pore.inv_Pc'])
        except:
            raise Exception('Orindary percolation has not been run!')
        tind = self._net.throats()
        conns = self._net.find_connected_pores(tind)
        for inv_val in inv_points[0:-1]:
            # Find clusters of defender pores
            Pinvaded = self['pore.inv_Pc'] <= inv_val
            Cstate = sp.sum(Pinvaded[conns], axis=1)
            Tinvaded = self['throat.inv_Pc'] <= inv_val
            # 0 = all open, 1=1 pore filled,
            # 2=2 pores filled 3=2 pores + 1 throat filled
            Cstate = Cstate + Tinvaded
            clusters = self._net.find_clusters(Cstate == 0)
            # Clean up clusters (invaded = -1, defended >=0)
            clusters = clusters * (~Pinvaded) - (Pinvaded)
            # Identify clusters connected to outlet sites
            out_clusters = sp.unique(clusters[p_outlets])
            trapped_pores = ~sp.in1d(clusters, out_clusters)
            trapped_pores[Pinvaded] = False
            if sum(trapped_pores) > 0:
                inds = (self['pore.trapped'] == 0) * trapped_pores
                self['pore.trapped'][inds] = inv_val
                trapped_throats = self._net.find_neighbor_throats(trapped_pores)
                trapped_throat_array = sp.asarray([False] * len(Cstate))
                trapped_throat_array[trapped_throats] = True
                inds = (self['throat.trapped'] == 0) * trapped_throat_array
                self['throat.trapped'][inds] = inv_val
                inds = (self['throat.trapped'] == 0) * (Cstate == 2)
                self['throat.trapped'][inds] = inv_val
        self['pore.trapped'][self['pore.trapped'] > 0] = 0
        self['throat.trapped'][self['throat.trapped'] > 0] = 0
        self['pore.inv_Pc'][self['pore.trapped'] > 0] = 0
        self['throat.inv_Pc'][self['throat.trapped'] > 0] = 0

    def return_results(self, Pc=0, seq=None, sat=None,
                       occupancy='occupancy_wp'):
        r"""
        Updates the occupancy status of invading and defending phases
        as determined by the OP algorithm
        """
        p_inv = self['pore.inv_Pc']
        p_inv[p_inv==sp.inf] = -1   # to be transferred to vtk
        self._phase_inv['pore.inv_Pc']=p_inv
        t_inv = self['throat.inv_Pc']
        t_inv[t_inv==sp.inf] = -1
        self._phase_inv['throat.inv_Pc']=t_inv
        #Apply invasion sequence values (to correspond with IP algorithm)
        p_seq = self['pore.inv_seq']
        self._phase_inv['pore.inv_seq']=p_seq
        t_seq = self['throat.inv_seq']
        self._phase_inv['throat.inv_seq']=t_seq
        #Apply saturation to pores and throats
        self._phase_inv['pore.inv_sat']=self['pore.inv_sat']
        self._phase_inv['throat.inv_sat']=self['throat.inv_sat']

        if(sat is not None):
            p_inv = self['pore.inv_sat'] <= sat
            t_inv = self['throat.inv_sat'] <= sat
            # Apply occupancy to invading phase
            temp = sp.array(p_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['pore.'+occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.'+occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.'+occupancy] = temp
                temp = sp.array(~t_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.'+occupancy] = temp
        elif(seq is not None):
            # definition: sat_i = seq_i/num_seq
            p_seq = self['pore.inv_seq'] <= seq
            t_seq = self['throat.inv_seq'] <= seq
            # Apply occupancy to invading phase
            temp = sp.array(p_seq, dtype=sp.float_, ndmin=1)
            self._phase_inv['pore.'+occupancy] = temp
            temp = sp.array(t_seq, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.'+occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_seq, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.'+occupancy] = temp
                temp = sp.array(~t_seq, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.'+occupancy] = temp
        else:
            p_inv = self['pore.inv_Pc'] <= Pc
            t_inv = self['throat.inv_Pc'] <= Pc
            # Apply occupancy to invading phase
            temp = sp.array(p_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['pore.'+occupancy] = temp
            temp = sp.array(t_inv, dtype=sp.float_, ndmin=1)
            self._phase_inv['throat.'+occupancy] = temp
            # Apply occupancy to defending phase
            if self._phase_def is not None:
                temp = sp.array(~p_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['pore.'+occupancy] = temp
                temp = sp.array(~t_inv, dtype=sp.float_, ndmin=1)
                self._phase_def['throat.'+occupancy] = temp

    def plot_moisture_curve(self,
                            pore_volume='volume',
                            throat_volume='volume',
                            pore_label='all',
                            throat_label='all'):
          r"""
          Plot drainage capillary pressure curve
          """
          try:
            PcPoints = sp.unique(self.get_data(prop='inv_Pc',pores='all'))
          except:
            raise Exception('Cannot print drainage curve: ordinary percolation simulation has not been run')
          pores=self._net.pores(labels=pore_label)
          throats = self._net.throats(labels=throat_label)
          Snwp_t = sp.zeros_like(PcPoints)
          Snwp_p = sp.zeros_like(PcPoints)
          Pvol = self._net['pore.'+pore_volume]
          Tvol = self._net['throat.'+throat_volume]
          Pvol_tot = sum(Pvol)
          Tvol_tot = sum(Tvol)
          for i in range(0,sp.size(PcPoints)):
              Pc = PcPoints[i]
              Snwp_p[i] = sum(Pvol[self._p_inv[pores]<=Pc])/Pvol_tot
              Snwp_t[i] = sum(Tvol[self._t_inv[throats]<=Pc])/Tvol_tot
          plt.plot(PcPoints,Snwp_p,'r.-')
          plt.plot(PcPoints,Snwp_t,'b.-')
          plt.xlim(xmin=0)
          plt.show()

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)
