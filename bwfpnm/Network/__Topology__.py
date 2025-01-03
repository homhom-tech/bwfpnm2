import scipy as sp
#from OpenPNM.Network.__GenericNetwork__ import GenericNetwork
from . import GenericNetwork
from bwfpnm.Utilities import transformations as trans
import OpenPNM.Utilities.misc as misc
from OpenPNM.Base import logging as _logging
from OpenPNM.Base import Controller as _controller
logger = _logging.getLogger(__name__)
ctrl = _controller()
#ctrl = OpenPNM.Base.Controller()


class Topology(GenericNetwork):
    '''Create topological network and assign labels to its pores.
    Labels assigned are internal, front, back, left, right, top, and bottom.

    example:
    import OpenPNM
    pn = OpenPNM.Network.Topoloy(name='net', coords=coords, conns=conns,
                                 isInlet=isInlet, isOulet=isOutlet)
    print(pn)'''
    def __init__(self, coords=[], conns=[], isInlet=False, isOutlet=False,
                 macro_Lx=None, macro_Ly=None, macro_Lz=None, **kwargs):
        super(Topology, self).__init__(**kwargs)

        Np = sp.shape(coords)[0]
        Nt = sp.shape(conns)[0]
        self['pore.coords'] = sp.array(coords)
        self['throat.conns'] = sp.array(conns)
        self['pore.all'] = sp.ones((Np,), dtype=bool)
        self['throat.all'] = sp.ones((Nt,), dtype=bool)
#        self['prop.all'] = sp.ones((Np,), dtype=bool)
        self['pore.index']   = sp.arange(0,len(self['pore.coords']))
        self['pore.inlet'] = isInlet
        self['pore.outlet'] = isOutlet

        self._macro_Lx = macro_Lx
        self._macro_Ly = macro_Ly
        self._macro_Lz = macro_Lz
        try:
            self._create_boundary_labels(['front', 'top', 'right', 'internal'])
        except:
            self['pore.internal'] = ~sp.logical_or(self['pore.inlet'],
                self['pore.outlet'])
        self._create_boundary_throats()

    def _create_boundary_labels(self,
                                labels=['front', 'top', 'right', 'internal']):
        x, y, z = self['pore.coords'].T
        xin, xout = x[self['pore.inlet']]  , x[self['pore.outlet']]
        eps = sp.amax([xin.max()-xin.min(), xout.max()-xout.min()])

        # X-direction: front-back
        if ('front' in labels) or ('back' in labels):
            try:    # for stitching purposes
                self['pore.front']
                self['pore.back']  = x <= x.min() + eps#x-rad <= x.min() #
                self['pore.front'] = x >= x.max() - eps#x+rad >= x.max() #
            except:
                self['pore.front'] = self['pore.outlet']
                self['pore.back']  = self['pore.inlet']
            self['pore.internal'] = self['pore.front'] + self['pore.back']
        else:
            self.del_properties(['pore.front', 'pore.back'])
        # Y-direction: right-left
        if ('right' in labels) or ('left' in labels):
            self['pore.left']  = y <= y.min()+eps#y-rad <= y.min()   #
            self['pore.right'] = y >= y.max()-eps#y+rad >= y.max()   #
            self['pore.internal'] += self['pore.left'] + self['pore.right']
        else:
            self.del_properties(['pore.left', 'pore.right'])
        # Z-direction: top-bottom
        if ('top' in labels) or ('bottom' in labels):
            self['pore.bottom'] = z <= z.min()+eps#z-rad <= z.min()  #
            self['pore.top']    = z >= z.max()-eps#z+rad >= z.max()  #
            self['pore.internal'] += self['pore.top'] + self['pore.bottom']
        else:
            self.del_properties(['pore.top', 'pore.bottom'])

        self['pore.internal'] = ~self['pore.internal']
        if 'internal' not in labels:
            self.del_properties(['pore.internal'])

    def _create_boundary_throats(self):
        tpin = self.find_neighbor_throats(self['pore.inlet'])
        tpout = self.find_neighbor_throats(self['pore.outlet'])
        self['throat.inlet'] = self.tomask(throats=tpin)
        self['throat.outlet'] = self.tomask(throats=tpout)

    def change_inoutlet(self, direction, update_bc=True):
        r'''
        direction: either x, y, or z.

        A GEOMETRY object is required if update_bc=True!
        '''
        # save the original inlet & outlet
        self['pore.inlet_ori'] = self['pore.inlet']
        self['pore.outlet_ori'] = self['pore.outlet']
        if update_bc:
            # redefined boundary pores
            self._geometries[0].update_boundary_labels()

        labels = {'x': ('back', 'front'), 'y': ('left', 'right'),
                  'z': ('bottom', 'top')}
        inlabel, outlabel = labels[direction]
        self._set_inoutlet(inlabel, outlabel, delete=False)

    def _set_inoutlet(self, inlabel, outlabel, delete=True):
        # save the original inlet & outlet
        self['pore.inlet'] = self['pore.'+inlabel]
        self['pore.outlet'] = self['pore.'+outlabel]
        self._create_boundary_throats()
        if delete:
            self.del_properties(['pore.'+inlabel, 'pore.'+outlabel])

    def _modify_bc_labels(self, Np, Nt, stitch_label='front', **kwargs):
        if stitch_label == 'inlet':
            stitch_label = 'back'
        elif stitch_label == 'outlet':
            stitch_label = 'front'

        labels_pair = {'front': 'back', 'back': 'front',
                       'top': 'bottom', 'bottom': 'top',
                       'right': 'left', 'left': 'right'}
        for key in labels_pair.keys():
            self['pore.'+key][Np:] = self['pore.'+key][:Np]

        self['pore.'+stitch_label][:Np] = False
        self['pore.'+stitch_label][Np:] = self['pore.'+labels_pair[stitch_label]][Np:]
        self['pore.'+labels_pair[stitch_label]][Np:] = False

        self['pore.inlet'] = self['pore.back']
        self['pore.outlet']  = self['pore.front']
        self._create_boundary_throats()

    def add_inoutlet(self, bc_tlength):
        r''' This method uses ``clone`` to clone the surface pores (labeled 'left',
        'right', etc), then shifts them to the periphery of the domain, and
        gives them the label 'right_face', 'left_face', etc.
        '''
        x, y, z = self['pore.coords'].T
        inletcoord = self['pore.coords'][self.pores('inlet')]
        xmean, ymean, zmean = sp.mean(inletcoord.T, axis=1)
        xyzdiff = sp.array([sum(abs(x-xmean)), sum(abs(y-ymean)), sum(abs(z-zmean))])
        xyzmin = xyzdiff - min(xyzdiff)

        # clone inlet & outlet pores
        for label in ['inlet', 'outlet']:
            ps = self.pores(label)
            self.clone(pores=ps, apply_label=[label+'_boundary', 'boundary'])

        self['throat.internal'] = self['throat.all'] - self['throat.boundary']
        # translate positions of cloned pores using bc throat length
        # bc_tlength is negative if connected to inlet
        bc_ind = self['pore.inlet'] | self['pore.outlet']
        bc_coords = self['pore.coords'][bc_ind].T
#        bc_tlength = bc_throat['lengthtotal']  # the same order with coords_bc provided sorted node & link.dat
        if xyzmin[0] < 1e-20:
            bc_coords[0] = bc_coords[0] + bc_tlength
        elif xyzmin[1] < 1e-20:
            bc_coords[1] = bc_coords[1] + bc_tlength
        else:
            bc_coords[2] = bc_coords[2] + bc_tlength
        self['pore.coords'][bc_ind] = bc_coords.T

    def domain_length(self, direction='x'):
        length = {'x': self._macro_Lx, 'y': self._macro_Ly, 'z': self._macro_Lz}
        return length[direction.lower()]

    def domain_area(self, face='xy'):
        length = {'x': self._macro_Lx, 'y': self._macro_Ly, 'z': self._macro_Lz}
        area = length[face[0].lower()]*length[face[1].lower()]
        return area

    def domain_bulk_volume(self):
        r'''Multiplying the dimensions: dx*dy*dz.
        '''
        vol = self._macro_Lx * self._macro_Ly * self._macro_Lz
        return vol

    def domain_pore_volume(self, pores=None, throats=None):
        r'''Adding all pore and throat volumes
        '''
        geo = self._geometries[0]
        try:
            pvol = sp.sum(geo['pore.volume'][pores])
        except:
            raise Exception("geo['pore.volume'] data is not available")
            pvol = 0
        try:
            tvol = sp.sum(geo['throat.volume'][throats])
        except:
            raise Exception("geo['throat.volume'] data is not available")
            tvol = 0
        return pvol+tvol

    def porosity(self, pores=None, throats=None):
        r'''Vpores/Vmaterials
        '''
        Vpore = self.domain_pore_volume(pores, throats)
        Vmat = self.domain_bulk_volume()
        return Vpore/Vmat

    def trim_geom_data(self, geoinput, trim_pores):
        r'''Modify geoinput with deleted trim_pores

        Arguments:
        ----------
        geoinput
        trim_pores
        '''
        trim_pores = sp.array(trim_pores,ndmin=1)
        Pkeep = sp.ones((self.num_pores(),),dtype=bool)
        Pkeep[trim_pores] = False
        Tkeep = sp.ones((self.num_throats(),),dtype=bool)
        Ts = self.find_neighbor_throats(trim_pores)
        if len(Ts)>0:
            Tkeep[Ts] = False
        keys = list(geoinput.keys())
        Tkeep2 = sp.array(list(zip(Tkeep, Tkeep)))
        for key in keys:
            if key[0] == 'p':
                geoinput[key] = geoinput[key][Pkeep]
            elif key[0] == 't':
                try:
                    geoinput[key] = geoinput[key][Tkeep]
                except:
                    geoinput[key] = geoinput[key][Tkeep2]

    def scaling(self, scalearray=[1,1,1], replace=False):
        r'''Scaling the pore coordinates and the global geometry
        by a scaling matrix [sx, sy, sz]. However,
        NOTE that only isotropic scaling is implemented for geometric properties!
        ==> scale = a constant

        Arguments:
        ----------
        scalearray     : an array of 3 scaling components for 3 directions x, y, and z: [sx, sy, sz].
                         This array is stored in self._scalearray.
        replace         : Boolean. True -> replace the network properties with the scaled ones. False -> return (coords, [Lx, Ly, Lz])
        '''
        if sp.size(scalearray) == 1:
            scalearray = sp.array([scalearray]*3)
        elif sp.size(scalearray) == 2:
            scalearray = sp.array(scalearray+[1])
        else:
            scalearray = sp.array(scalearray)

        coords = self['pore.coords']*scalearray
        Lx = self._macro_Lx*scalearray[0]
        Ly = self._macro_Ly*scalearray[1]
        Lz = self._macro_Lz*scalearray[2]
        if replace:
            self['pore.coords'] = coords
            self._macro_Lx = Lx
            self._macro_Ly = Ly
            self._macro_Lz = Lz
            self._scalearray = scalearray
            return
        else:
            return (coords, [Lx, Ly, Lz])

    def cluster_types(self, mask=[], save=True):
        r'''Identifying cluster types: spanning, deadend, and isolated.
        If save=True: Create network properties: pore.span, pore.deadend,
        pore.isolated.

        Spanning clusters: pores and throats connected to both inlet and outlet
        Dead-end clusters: pores and throats that are not spanning nor isolated.
        Isolated clusters: pores and throats not connected to both inlet and outlet
        Surface clusters: pores and throats connected to either inlet or oulet.

        Return:     Array of (spanning, deadend, isolated) clusters.
                    Each cluster is an array of pores.
        -------
        '''
        # Check for separated clusters of pores
        disc_Cs = []
        temp = []
        if len(mask)==0:
            mask = self.tomask(throats=self.throats('all'))
        # Cs: an Np list of cluster numbers, started from 0 to Nc-1
        Cs = self.find_clusters(mask)
        for i in sp.unique(Cs):
            temp.append(sp.where(Cs == i)[0]) # Nc list of pore id
        b = sp.array([len(item) for item in temp]) # Nc list of cluster sizes
        c = sp.argsort(b)[::-1] # Nc sorted list of arg from biggest cluster
        for i in range(0, len(c)):
            disc_Cs.append(temp[c[i]]) # Nc list: sorted temp w.r.t. size

        self['pore.cluster_id'] = Cs
        self['throat.cluster_id'] = Cs[self['throat.conns']].mean(axis=1)
#        allpores = self.pores()

        clusters = sp.array(disc_Cs)
        inlet = self['pore.inlet']
        outlet = self['pore.outlet']
        isbc = []
        for cluster in clusters:
            isinlet = sp.any(inlet[cluster])
            isoutlet = sp.any(outlet[cluster])
            isbc.append((isinlet, isoutlet))
        isbc = sp.array(isbc)
        spanning_clusters = sp.all(isbc, axis=1)  #isinlet=True & isoutlet=True
        isolated_clusters = sp.all(~isbc, axis=1) #isinlet=0 & isoutlet=0
        deadend_clusters = ~(spanning_clusters+isolated_clusters) # = surface without spanning
        surface_clusters = sp.any(isbc, axis=1) #isinlet=True or isoutlet=True (spanning included)
        inlet_clusters = isbc[:,0] * deadend_clusters # isinlet=1 & isoutlet=0
        outlet_clusters = isbc[:,1] * deadend_clusters # isinlet=0 & isoutlet=1

        # pores belong to the following clusters:
        span = clusters[spanning_clusters]
        dead = clusters[deadend_clusters]
        isolated = clusters[isolated_clusters]
        surface = clusters[surface_clusters]
        inletc = clusters[inlet_clusters]
        outletc = clusters[outlet_clusters]

        try:
            span = sp.hstack(span)
            spant = self.find_neighbor_throats(pores=span)
        except:
            span, spant = [], []
        try:
            dead = sp.hstack(dead)
            deadt = self.find_neighbor_throats(pores=dead)
        except:
            dead, deadt = [], []
        try:
            pisol = sp.hstack(isolated)
            tisol = self.find_neighbor_throats(pores=pisol)
        except:
            pisol, tisol = [], []

        if save:
            try:
                pinletc = sp.hstack(inletc)
                tinletc = self.find_neighbor_throats(pores=pinletc)
            except:
                pinletc, tinletc = [], []
            try:
                poutletc = sp.hstack(outletc)
                toutletc = self.find_neighbor_throats(pores=poutletc)
            except:
                poutletc, toutletc = [], []
            self['pore.span'] = self.tomask(pores=span)
            self['pore.surface'] = self.tomask(pores=dead) #renamed to 'surface'
            self['pore.deadend'] = self['pore.surface']
            self['pore.isolated'] = self.tomask(pores=pisol)
#            self['pore.surface'] = self.tomask(pores=clusters[surface][0])
            self['pore.inlet_clusters'] = self.tomask(pores=pinletc)
            self['pore.outlet_clusters'] = self.tomask(pores=poutletc)
            # related TRHOATS:
#            surfacet = self.find_neighbor_throats(pores=self.pores('surface'))
            self['throat.span'] = self.tomask(throats=spant)
            self['throat.surface'] = self.tomask(throats=deadt)
            self['throat.deadend'] = self['throat.surface']
            self['throat.isolated'] = self.tomask(throats=tisol)
#            self['throat.surface'] = self.tomask(throats=surfacet)
            self['throat.inlet_clusters'] = self.tomask(throats=tinletc)
            self['throat.outlet_clusters'] = self.tomask(throats=toutletc)

        return (span, dead, pisol)

    def span_existence(self, mask=[], return_isolated=True):
        r'''Identifying whether a spanning cluster is formed,
        if so return the cluster.

        Arguments
        ---------
        mask: array_like, boolean
            A list of active nodes.  This method will automatically search
            for clusters based on site or bond connectivity depending on
            wheather the received mask is Np or Nt long.

        Return:
        If spanning cluster exists return: (True, spanning cluster)
        Else: (False, one biggest cluster)
        -------
        '''
        # Check for separated clusters of pores
        disc_Cs = []
        temp = []
        Cs = self.find_clusters(mask=mask)
        for i in sp.unique(Cs):
            temp.append(sp.where(Cs == i)[0])
        b = sp.array([len(item) for item in temp])
        c = sp.argsort(b)[::-1]
        for i in range(0, len(c)):
            disc_Cs.append(temp[c[i]])

        clusters = sp.array(disc_Cs)
        inlet = self['pore.inlet']
        outlet = self['pore.outlet']
        isbc = []
        for cluster in clusters:
            isinlet = sp.any(inlet[cluster])
            isoutlet = sp.any(outlet[cluster])
            isbc.append((isinlet, isoutlet))
        isbc = sp.array(isbc)
        spanning_clusters = sp.all(isbc, axis=1)
        isolated_clusters = sp.all(~isbc, axis=1)
#        deadend_clusters = ~(spanning_clusters+isolated_clusters)

        span = sp.where(spanning_clusters)[0]
#        dead = sp.where(deadend_clusters)[0]
        isolated = sp.where(isolated_clusters)[0]

        existence = bool(span.size)
        if existence:
            if return_isolated:
                return (existence, clusters[span], clusters[isolated])
            else:
                return (existence, clusters[span])
        else:
            return (existence, clusters[0])

    def stitch(self, donor, P_network, P_donor, method='nearest',
               len_max=sp.inf, label_suffix=''):
        r'''
        ---This is a customized version of OpenPNM.Utilities.topology.stich()--
        -- modification: stitch only point to point (not a combination)

        Stitches a second a network to the current network.

        Parameters
        ----------
        networK : OpenPNM Network Object
            The Network that will to which to donor Network will be attached

        donor : OpenPNM Network Object
            The Network to stitch on to the current Network

        P_network : array_like
            The pores on the current Network

        P_donor : array_like
            The pores on the donor Network

        label_suffix : string or None
            Some text to append to each label in the donor Network before
            inserting them into the recipient.  The default is to append no
            text, but a common option would be to append the donor Network's
            name. To insert none of the donor labels, use None.

        len_max : float
            Set a length limit on length of new throats

        method : string (default = 'delaunay')
            The method to use when making pore to pore connections. Options are:

            - 'delaunay' : Use a Delaunay tessellation (not implemented)
            - 'nearest' : Connects each pore on the receptor network to its nearest
                          pore on the donor network

        Notes
        -----
        Before stitching it is necessary to translate the pore coordinates of
        one of the Networks so that it is positioned correctly relative to the
        other.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = OpenPNM.Network.TestNet()
        >>> [pn.Np, pn.Nt]
        [125, 300]
        >>> [pn2.Np, pn2.Nt]
        [125, 300]
        >>> pn2['pore.coords'][:, 2] += 5.0
        >>> pn.stitch(donor=pn2, P_network=pn.pores('top'),
        ...           P_donor=pn2.pores('bottom'), method='nearest', len_max=1.0)
        >>> [pn.Np, pn.Nt]
        [250, 625]

        '''
        # Ensure Networks have no associated objects yet
        if (len(self._simulation()) > 1) or (len(donor._simulation()) > 1):
            raise Exception('Cannot stitch a Network with active sibling objects')
        # Get the initial number of pores and throats
        N_init = {}
        N_init['pore'] = self.Np
        N_init['throat'] = self.Nt
        if method == 'nearest':
            P1 = P_network
            P2 = P_donor + N_init['pore']  # Increment pores on donor
            C1 = self['pore.coords'][P_network]
            C2 = donor['pore.coords'][P_donor]
            D = sp.linalg.norm(C1-C2, axis=1)   # MODIFICATON here!
            P_ind = sp.where(D <= len_max)[0]
            conns = sp.vstack((P1[P_ind], P2[P_ind])).T
        else:
            raise RuntimeError('<{}> method not supported'.format(method))

        # Enter donor's pores into the Network
        self.extend(pore_coords=donor['pore.coords'])
        self['pore.index'] = sp.arange(self.Np)

        # Enter donor's throats into the Network
        self.extend(throat_conns=donor['throat.conns'] + N_init['pore'])

        # Add all donor labels to recipient network
        if label_suffix is not None:
            if label_suffix != '':
                label_suffix = '_'+label_suffix
            for label in donor.labels():
                element = label.split('.')[0]
                locations = sp.where(self._get_indices(element) >=
                                      N_init[element])[0]
                try:
                    self[label + label_suffix]
                except:
                    self[label + label_suffix] = False
                self[label+label_suffix][locations] = donor[label]

        # Add the new stitch throats to the Network
        self.extend(throat_conns=conns, labels='stitched')

        # Remove donor from Controller, if present
        # This check allows for the reuse of a donor Network multiple times
        if donor in ctrl.values():
            ctrl.purge_object(donor)

    def clone_network(self, stitch_label='outlet', apply_label=['clone'],
                      clone_geo=True, **kwargs):
        r'''
        Clones the whole network (pores & throats) and adds them to the network
        by mirroring to the stitch_label pores.

        Parameters
        ----------
        apply_labels : string, or list of strings
            The labels to apply to the clones, default is 'clone'
        stitch_label : string
            pore label to be stitched
        '''
        if (self._geometries != []):
            logger.warning('Network has active Geometries, new pores must be \
                            assigned a Geometry')
        if (self._phases != []):
            raise Exception('Network has active Phases, cannot proceed')

        self._stitch_label = stitch_label

        logger.debug('Cloning pores')
        apply_label = list(apply_label)

        # Clone pores
        Np = self.Np    # copy the original size
        Nt = self.Nt
        pcoords = self['pore.coords']

        # reflect clone coordinates
        M = self._create_reflection_matrix(label=stitch_label)
        pclone = self._reflect_coords(M, pcoords)

        # Enter clone's pores into the Network
        self.extend(pore_coords=pclone)
        self['pore.index'] = sp.arange(2*Np)
        self['pore.original'] = self['pore.all']
        self['pore.original'][Np:] = False

        # Add/initialize clone labels to network
        for item in apply_label:
            if 'pore.' + item not in self.keys():
                self['pore.'+item] = False
            if 'throat.' + item not in self.keys():
                self['throat.'+item] = False

        # Enter clone's throats into the Network
        tclone = self['throat.conns'] + Np
        self.extend(throat_conns=tclone)
        self['throat.original'] = self['throat.all']
        self['throat.original'][Nt:] = False

        # Apply provided labels to cloned pores & throats
        for item in apply_label:
            self['pore.'+item][self.pores('all') >= Np] = True
            self['throat.'+item][self.throats('all') >= Nt] = True

        # stitch parent & clone networks
        pores1 = self.pores(stitch_label)
        pores2 = pores1 + Np
        conns = self._stitch_pores(pores1, pores2,
                                   pcoords[pores1], pclone[pores1])

        # Add the new stitch throats to the Network
        self.extend(throat_conns=conns, labels='stitched')
        t_stitched_new = self.throats('stitched')[-conns.shape[0]:]

        # LABEL ADJUSMENT!
        self._modify_bc_labels(Np, Nt, stitch_label)

        # Any existing adjacency and incidence matrices will be invalid
        self._update_network()

        # Clone geometry as well
        if clone_geo:
            self._geometries[0].clone_geometry(t_stitched_new)

    def _stitch_pores(self, pores1, pores2, coords1, coords2, len_max=sp.inf):
        D = sp.linalg.norm(coords1-coords2, axis=1)   # MODIFICATON here!
        P_ind = sp.where(D <= len_max)[0]
        conns = sp.vstack((pores1[P_ind], pores2[P_ind])).T
        return conns

    def _create_reflection_matrix(self, label='outlet', distance=None):
        '''Determine reflection plane (point, normal), then the reflection matrix M
        '''
        pplane = self.pores(label)
        pcoords = self['pore.coords']

        axis_pos = sp.argmin(sp.std(pcoords[pplane], axis=0))
        mean_pos = sp.mean(pcoords[:, axis_pos])
        if distance is None:
            # set the outest pore to plane distance = a*lengthtotal.mean
            distance = sp.mean(sp.absolute(sp.diff(pcoords[:, axis_pos])))/10

        if sp.all(pcoords[pplane][:, axis_pos] > mean_pos):
            ref_pos = sp.argmax(pcoords[pplane][:, axis_pos])
        elif sp.all(pcoords[pplane][:, axis_pos] < mean_pos):
            ref_pos = sp.argmin(pcoords[pplane][:, axis_pos])
            distance *= -1
        else:
            raise Exception('Reference point for reflection cannot be determined')

        point = pcoords[pplane][ref_pos]
        point[axis_pos] += distance
        normal = sp.zeros(3)
        normal[axis_pos] = 1
        M = trans.reflection_matrix(point, normal)
        return M

    def _reflect_coords(self, M, coords):
        # add the forth coordinate equals 1
        pcoord2 = sp.ones((coords.shape[0], 4))
        pcoord2[:, :3] = coords
        # reflect them
        pcoord2_new = sp.dot(M, pcoord2.T)
        coords = pcoord2_new[:3, :].T
        return coords

    def del_properties(self, props=[]):
        for prop in props:
            try:
                del self[prop]
            except:
                pass

    def network_properties(self, netname='Berea', resolution=5.345):
        geo = self._geometries[0]
        L = [self._macro_Lx, self._macro_Ly, self._macro_Lz]
        porosity = (geo['pore.volume'].sum()+geo['throat.volume'].sum())/(sp.cumprod(L)[-1])
        geo.count_shape()
        Ncirc = geo['pore.circular'].sum()+geo['throat.circular'].sum()
        Nsqu = geo['pore.square'].sum()+geo['throat.square'].sum()
        Ntri = geo['pore.triangular'].sum()+geo['throat.triangular'].sum()
        # Average coordination number
        z = self.num_neighbors(pores=self.pores())
        z_avg = z.sum()/z.size

#        try:
#            self['pore.span']
#        except:
        span, surface, isolated = self.cluster_types(save=True)

        props = {'Resolution': resolution,
                 'Porosity': porosity,
                 'Side length': L[0],
                 'N-pores': self.Np,
                 'N-throats': self.Nt,
                 'N-clusters': sp.unique(self['pore.cluster_id']).size,
                 'N-pore-inlets': self.pores('inlet').size,
                 'N-pore-outlets': self.pores('outlet').size,
                 'N-circular': Ncirc,
                 'N-square': Nsqu,
                 'N-triangular': Ntri,
                 'Pore radius-avg': geo['pore.diameter'].mean()/2,
                 'Pore radius-min': geo['pore.diameter'].min()/2,
                 'Pore radius-max': geo['pore.diameter'].max()/2,
                 'Throat radius-avg': geo['throat.diameter'].mean()/2,
                 'Throat radius-min': geo['throat.diameter'].min()/2,
                 'Throat radius-max': geo['throat.diameter'].max()/2,
                 'Pore length-avg': geo['throat.porelengths'].mean(),
                 'Pore length-min': geo['throat.porelengths'].min(),
                 'Pore length-max': geo['throat.porelengths'].max(),
                 'Throat length-avg': geo['throat.length'].mean(),
                 'Throat length-min': geo['throat.length'].min(),
                 'Throat length-max': geo['throat.length'].max(),
                 'G-pore avg': geo['pore.shapefactor'].mean(),
                 'G-pore min': geo['pore.shapefactor'].min(),
                 'G-pore max': geo['pore.shapefactor'].max(),
                 'G-throat avg': geo['throat.shapefactor'].mean(),
                 'G-throat min': geo['throat.shapefactor'].min(),
                 'G-throat max': geo['throat.shapefactor'].max(),
                 'Coord number-avg': z_avg,
                 'Coord number-min': z.min(),
                 'Coord number-max': z.max(),
                 'N-spanning pores': self['pore.span'].sum(),
                 'N-surface pores': self['pore.surface'].sum(),
                 'N-isolated pores': self['pore.isolated'].sum(),
                 'N-spanning throats': self['throat.span'].sum(),
                 'N-surface throats': self['throat.surface'].sum(),
                 'N-isolated throats': self['throat.isolated'].sum()}
        try:
            props.update({'N-fine pores': self['pore.fine'].sum()})
        except:
            props.update({'N-fine pores': 'undefined'})

        try:
            props.update({'N-fine throats': self['throat.fine'].sum()})
        except:
            props.update({'N-fine throats': 'undefined'})

        hline = '-' * 60
        lines = [hline]
        lines.append('Network properties: \t ' + netname)
        lines.append(hline)
        lines.append("{0:<5s} {1:<35s} {2:<10s}".format('#', 'Properties', 'Values'))
        lines.append(hline)
        for i, item in enumerate(props.keys()):
            prop = item
            if len(prop) > 35:
                prop = prop[0:32] + '...'
            try:
                lines.append("{0:<5d} {1:<35s} {2:>5g}".format(i + 1,
                         prop, props[item]))
            except:
                pass
        lines.append(hline)
        print('\n'.join(lines))

    def fluid_properties(self):
        hline = '-' * 60
        lines = [hline]
        for fluid in self._phases:
            lines.append('Fluid properties: \t {}'.format(fluid.name))
            lines.append(hline)
            lines.append("{0:<5s} {1:<35s} {2:<10s}".format('#', 'Properties', 'Values'))
            lines.append(hline)
            for i, item in enumerate(fluid.props()):
                prop = item
                if len(prop) > 35:
                    prop = prop[0:32] + '...'
                try:
                    lines.append("{0:<5d} {1:<35s} {2:>5g}".format(i + 1,
                             prop, fluid[item][0]))
                except:
                    pass
            lines.append(hline)
        print('\n'.join(lines))

    def iscoplanar(self, plane):
        r'''
        plane: indices of the pores
        '''
        pcoords = self['pore.coords'][plane]
        return misc.iscoplanar(pcoords)

    def fine_elements(self, Dmin, case=2, save=True):
        r'''
        Dmin: Threshold diameter below which are fine elements
        Case:
            1:  Only fine pores
            2:  Only fine throats,
                but preserve inlet/outlet throats
            3:  Fine pores & throats,
                but preserve fine pores connecting 2 coarse pores

        Return:     (fine pores, fine throats) in boolean lists
        '''
        Pfine = self['pore.diameter'] < Dmin
        Tfine = self['throat.diameter'] < Dmin
        if case == 1:   # Fine pores only
            # Pfine = Pfine
            Tfine = False
        elif case == 2: # Fine throats only (resulting isolated fine pores)
            # but, this creates problem in permeability calc, due to
            # isolated/single inlet/outlet pores.
            # So, keep the fine inlet/outlet throats!
            Pfine = False
            Tidsave = sp.r_[self.throats('inlet'), self.throats('outlet')]
            Tfine[Tidsave] = False
        elif case == 3: # Fine throats & fine pores
            Tcoarse = ~Tfine
            Pidsave = []
            # preserved Pfine if Tcoarse-Pfine-Tcoarse
            Pfineid = self.toindices(Pfine)
            temp = self.find_neighbor_throats(pores=Pfineid, flatten=False)
            for i, throat in list(zip(Pfineid,temp)):
                try:
                    if Tcoarse[throat].sum() >= 2:
                        Pidsave.append(i)
                except:
                    pass
            Pfine[Pidsave] = False

        if save:
            self['pore.fine'] = Pfine
            self['throat.fine'] = Tfine

        return (Pfine, Tfine)

    def compare_objs(obj1, obj2):
        r'''
        Compare the properties and their values of two objects.
        '''
        result = True
        keys1 = list(obj1.keys())
        keys2 = list(obj2.keys())
        keys1.sort()
        keys2.sort()
        try:
            assert sp.all(keys1==keys2)
        except:
            result = False
            print('The property names are not the same.')
        for item in keys1:
            try:
                assert sp.all(obj1[item]==obj2[item])
            except:
                result = False
                print('The values of {} are not the same'.format(item))
        return result

    def cut_network(self, ratio=0.5, center=[0.5, 0.5, 0.5], shift=[0, 0, 0],
                    flowdir='all', replace=True, suffix=''):
        r'''
        Cut the original network dimension based on the given ratio.
        Boundary pores and throats are added by modifying the original throats
        passing through the new boundary planes in the desired flow direction.

        Args:
        =================
        ratio: float (0, 1]
        center: [rx, ry, rz] determines the location of the new domain:
                r in [0, 1]: r=0   --> new_min = old_min (left)\n
                             r=0.5 --> new_mean = old_mean (center)\n
                             r=1   --> new_max = old_max (right)\n
        shift: [tx, ty, ty] --> shift the new domain (not impelemented)
                t in [-1, 1]
        flowdir: 'all' (default = ['x', 'y', 'z']), 'x', 'y', 'z'
        '''
        if ratio == 1:
            return

        pcoords = self['pore.coords']
        center = sp.array(center)
        shift = sp.array(shift)

        #%% Determine the new domain & pores
        Lnet = sp.array([self._macro_Lx, self._macro_Ly, self._macro_Lz])
        Lnet2 = Lnet*ratio
        # pmin2 and pmax2 determine the new domain/box
        Lspare = Lnet - Lnet2
        pmin2 = sp.amin(pcoords, axis=0) + Lspare*center
        pmax2 = pmin2 + Lnet2
        if sp.any(shift):
            pass # can be represented by the center alone

        xcond1 = pcoords[:, 0] > pmin2[0] # inplane pores are excluded here
        xcond2 = pcoords[:, 0] < pmax2[0] # so that those are identified as offbc pores
        ycond1 = pcoords[:, 1] > pmin2[1]
        ycond2 = pcoords[:, 1] < pmax2[1]
        zcond1 = pcoords[:, 2] > pmin2[2]
        zcond2 = pcoords[:, 2] < pmax2[2]
        # determine the new domain
        xcond = sp.logical_and(xcond1, xcond2)
        ycond = sp.logical_and(ycond1, ycond2)
        zcond = sp.logical_and(zcond1, zcond2)
        # determine the pores inside the new domain
        pin = sp.logical_and(xcond, ycond) # Pores INside the domain
        pin = sp.logical_and(pin, zcond) # True if it is in the new domain
        pinid = sp.where(pin)[0]

        #%% Identify the boundary pores & throats
        tconns = self['throat.conns']
        Nt = self.num_throats()
        # BC throats
        cond1 = sp.in1d(tconns[:, 0], pinid)
        cond2 = sp.in1d(tconns[:, 1], pinid)
        tbcnot = sp.logical_and(cond1, cond2) #inbc-to-inbc throats excluded
        cond1 = sp.logical_xor(cond1, tbcnot)
        cond2 = sp.logical_xor(cond2, tbcnot)
        idcond1 = sp.where(cond1)[0]
        idcond2 = sp.where(cond2)[0]
        tbc = sp.logical_xor(cond1, cond2)  # mask of boundary throats
        tbcid = sp.where(tbc)[0]            # indices of boundary throats
        # non-unique out and inside bcpores to conserve the connection
        idcond12 = sp.r_[idcond1, idcond2]
        argcond = sp.argsort(idcond12) # to match the order of tbcid
        idcond12 = idcond12[argcond]
        # indices of boundary pores: in & out of the domain
        # IN: pinbc_id - tbcid - poffbc_id :OUT
        pinbc_id = sp.r_[tconns[cond1][:,0], tconns[cond2][:,1]][argcond]
        poffbc_id = sp.r_[tconns[cond1][:,1], tconns[cond2][:,0]][argcond]

        #%% Force in-plane boundary pores - in-bcplane coordinates of pnewid
        # coordinates (x,y,z) of boundary pores (forced to be in the in/outplanes)
        cen = self.center_point(pmin2, pmax2)
        planes = ['x', 'y', 'z']
        normals = [sp.array([1,0,0]), sp.array([0,1,0]), sp.array([0,0,1])]
        points = [[sp.array([pmin2[0], cen[1], cen[2]]), #a point in x-min plane
                   sp.array([pmax2[0], cen[1], cen[2]])], #a point in x-max plane
                  [sp.array([cen[0], pmin2[1], cen[2]]),
                   sp.array([cen[0], pmax2[1], cen[2]])],
                  [sp.array([cen[0], cen[1], pmin2[2]]),
                   sp.array([cen[0], cen[1], pmax2[2]])]]
        pcords = dict(zip(planes, points))
        normaxis = list(zip(normals, planes))

        # Calculate the bc pores coordinates (reposition the pnewid)
        p0 = pcoords[pinbc_id]
        p1 = pcoords[poffbc_id]
        p1new = sp.zeros_like(p0) # in-bcplane coordinates
        p0bc = sp.zeros_like(p0)
        p1bc = sp.zeros_like(p0)
#        pisect = {'x':{}, 'y':{}, 'z':{}}
        for normal, axis in normaxis: # define the axis
            for i, pplane in enumerate(pcords[axis]): # define the min-max side
                psect, inds = self._isect_line_plane(p0, p1, normal, pplane)
                p1new[inds] = psect
#                p0bc[inds], p1bc[inds] = self._dist_to_bc(p0[inds], p1[inds],
#                                                          psect, normal)

        #%% Analyse the cases
        geo = self._geometries[0]
        rad = geo['pore.diameter']/2
        plength = geo['throat.porelengths']
        Lpinbc = sp.r_[plength[idcond1][:,0], plength[idcond2][:,1]][argcond]
        Lpoffbc = sp.r_[plength[idcond1][:,1], plength[idcond2][:,0]][argcond] #offbc, not the new pores

#        Din = sp.linalg.norm(p0bc- p0, axis=1)
#        Doff = sp.linalg.norm(p1bc - p1, axis=1)
#        case1 = sp.logical_or(Din <= rad[pinbc_id], Din <= Lpinbc)
#        casetemp = sp.logical_or(Doff<=rad[poffbc_id], Doff<=Lpoffbc)
#        case2 = sp.logical_and(~case1, casetemp)
#        case3 = sp.logical_and(~case1, ~casetemp)
##        case4 = sp.zeros_like(case3, dtype=bool)
#        idcase3 = sp.where(case3)[0]
#        for i, ip in list(zip(idcase3, p0[case3])):
##            Dr = sp.linalg.norm(ip - p0[case3], axis=1)
##            case4 = Dr < sp.amax([rad, ])
#            Dr = sp.zeros_like(Din)
#            Dd = sp.zeros_like(Din)
#            case4 = sp.zeros_like(case3, dtype=bool)
#            case5 = sp.zeros_like(case3, dtype=bool)
#            for j, jp in list(zip(idcase3, p0[case3])):
#                if i != j:
#                    Dr[j] = sp.linalg.norm(ip - jp, axis=1)
#                    case4[j] = Dr[j] < sp.amax([rad[i], rad[j]])
#                    case5[j] = Dr[j] < rad[i] + rad[j]
#            # replace
#
##        for i in sp.arange(tbcid.size):

        #%% Clone the inbc pores as many as their throats, as the intersection pores
        # in the boundary planes.
        # This modifies the network sizes: Npores and Nthroats
        label = 'bc'
        endlabel = '_0p' + str(ratio).split('.')[1] + suffix
        if not replace:
            label += endlabel
        self.clone_pores(pores=pinbc_id, apply_label=[label], mode='parents')
        pnewid = sp.arange(pin.size, self.Np) #index of the clone pores
        tnewid = sp.arange(Nt, self.Nt)
        pin = sp.concatenate((pin, sp.ones(pinbc_id.size, dtype=bool)))
        labels = {'pbc': sp.zeros_like(pin, dtype=bool),
                  'tbc': sp.zeros(self.Nt, dtype=bool),
                  'pinlet': sp.zeros_like(pin, dtype=bool),
                  'poutlet': sp.zeros_like(pin, dtype=bool),
                  'tbcid': tbcid,
                  'pinbcid': pinbc_id,
                  'poffbcid': poffbc_id}
        labels['pbc'][pnewid] = True
        labels['tbc'][tnewid] = True

        self['pore.coords'][pnewid] = p1new  # coords of the clone pores (pnewid)
        #%% Modify the geometrical properties of the boundary pores & throats
        # volume-wise: the resized bcthroats = bcthroats + bcpores
        # size-wise: bcthroats = bcpores
        # Extend geometrical properties
        geo = self._geometries[0]
        # Adjust/extend 'all' labels to the new size
        del geo['pore.all'], geo['throat.all']
        try:
            del self['pore.geo'], self['throat.geo']
        except:
            pass
        geo['pore.all'] = sp.ones((self.Np,), dtype=bool)
        geo['throat.all'] = sp.ones((self.Nt,), dtype=bool)
        self['pore.geo'] = True
        self['throat.geo'] = True
        for prop in geo.props('pore'):
            geo[prop] = sp.concatenate((geo[prop], geo[prop][poffbc_id]), axis=0)
        for prop in geo.props('throat'):
            geo[prop] = sp.concatenate((geo[prop], geo[prop][tbcid]), axis=0)

        # replace: bcpore sizes = bcthroat sizes
        props = ['diameter', 'shapefactor', 'shapefactor_constant', 'area']
        for prop in props:
            geo['pore.'+prop][pnewid] = geo['throat.'+prop][tbcid]

        ## ---- start: check data consistency ------


        # debug: the old conduit lengths
        Lconduitold = plength[tbcid].sum(axis=1)+geo['throat.length'][tbcid]
        Lconduitold2 = geo['throat.length'][tbcid]+Lpinbc+Lpoffbc #check Lpinbc & Lpoffbc
        Lconduitold3 = geo['throat.lengthtotal'][tbcid] # check data consistency:
        dL2 = Lconduitold - Lconduitold2
        dL3 = Lconduitold - Lconduitold3
        ## ---- end: check data consistency ------

        # LENGTHS of the NEW pores (pnewid) and throats (tnewid)
        # -------------------------------------------------------
        Lconduit = sp.linalg.norm(p1 - p0, axis=1) # new p1 - p0
        Lptbc = Lconduit - Lpinbc  # Lpoffbc + Ltbc
        mask = Lptbc > 0   # accepted pnewid: pnewid outside the pinbc_id
        Lpbc = Lptbc[mask] * 1/10     # poffbc_id lengths
        Ltbc = Lptbc[mask] * 9/10     # tbc lengths
#        Lpinbc = Lpinbc[mask]

        # old lengths
        Ltbcold = geo['throat.length'][tbcid][mask]
#        Lpinbcold = Lpinbc
#        Lpbcold = Lpoffbc[mask]
#        Lconduitold = Ltbcold + Lpinbcold + Lpbcold
#        Lptbcold = Ltbcold + Lpbcold

#        # Check the length data
#        sp.all(Lconduitold - Lconduit[mask] > 0)
#        sp.all(Lptbcold - Lptbc[mask] > 0)
#        sp.all(Lpbcold - Lpbc > 0)
#        sp.all(Lpinbcold - Lpinbc == 0) # the same
#        sp.all(Ltbcold - Ltbc > 0)

        # ADJUST the volumes accordingly
        # bugs: ratio > 1 !
        tratio = Ltbc / Ltbcold #Lnew/Lold
        tvolume = geo['throat.volume'][tbcid[mask]] * tratio *9/10
        pvolume = geo['throat.volume'][tbcid[mask]] * tratio *1/10

        # Special CASE: Ltbcold < Ltbc
        ids = sp.where(Ltbcold < Ltbc)[0]
        Ltbc[ids] = Ltbcold[ids]
        Lpbc[ids] = Lptbc[mask][ids] - Ltbc[ids]
        tvolume[ids] = geo['throat.volume'][tbcid[mask]][ids]
        pvolume[ids] = geo['pore.volume'][poffbc_id[mask]][ids]

        geo['throat.lengthtotal'][tnewid[mask]] = Lconduit[mask]
        geo['throat.length'][tnewid[mask]] = Ltbc
        geo['throat.volume'][tnewid[mask]] = tvolume
        geo['pore.volume'][pnewid[mask]] = pvolume

        # Update the new in-plane bcpore lengths
        plengthnew = geo['throat.porelengths'][:,1]
        plengthnew[tnewid[mask]] = Lpbc
#        plengthnew[:,1][cond2[tbcid[mask]]] = Lpbc[cond2[tbcid[mask]]]
#        plengthnew[:,1][cond1[tbcid[mask]]] = Lpbc[cond1[tbcid[mask]]]
        geo['throat.porelengths'][:,1] = plengthnew

        #%% Filter the pnewid:
        # delete the poffbc if (pinbc intersects the boundary):
        # 1) Lconduit - length_pinbc < 0  or
        # 2) Lconduit - rad_pinbc < 0
        # then adjust the properties of the pinbc
        offmask1 = sp.where(~mask)[0] # 1) Lconduit - length_pinbc < 0
        pradin = geo['pore.diameter'][pinbc_id]/2
        mask2 = Lconduit - pradin > 0
        offmask2 = sp.where(~mask2)[0] # 2) Lconduit - rad_pinbc < 0
        # choose the smallest ratio and adjust the volume of pbcin
        offmask = sp.unique(sp.r_[offmask1, offmask2])
        ratiol = Lconduit[offmask]/Lpinbc[offmask]
        ratior = Lconduit[offmask]/pradin[offmask]
        prat = sp.amin([ratiol, ratior], axis=0)

        pnewok = sp.logical_and(mask, mask2)
        pin[pnewid[offmask]] = False   # delete new bcpores included in pbcin
#        pkeepid = sp.r_[pinid, pnewid[~offmask]]

        #%% Create labels
        # gfg
        self['pore.inBC' + endlabel] = False
        self['pore.inBC' + endlabel][pinbc_id[pnewok]] = True
        self['pore.offBC' + endlabel] = False
        self['pore.offBC' + endlabel][poffbc_id] = True
        self['pore.BC' + endlabel] = False
        self['pore.BC' + endlabel][pnewid[pnewok]] = True
        self['pore.BC' + endlabel][pinbc_id[~pnewok]] = True

        #%% Apply the boundary labels
        # Initiate empty labels
        elements = ['pore.', 'throat.']
        for elem in elements:
            for temp in planes:
                if not replace:
                    temp += '_0p' + str(ratio).split('.')[1]
                self[elem + 'inlet_'+temp] = False
                self[elem + 'outlet_'+temp] = False
        # update the labels
        pco = self['pore.coords']
        for i, axis in enumerate(planes):
            temp1 = sp.where(pco[:,i] == pmin2[i])[0]
            temp2 = sp.where(pco[:,i] == pmax2[i])[0]

            for side, temp in list(zip(['inlet', 'outlet'], [temp1, temp2])):
                prop = axis
                if not replace:
                    prop += '_0p' + str(ratio).split('.')[1]
                tempt = self.find_neighbor_throats(temp)
                self['pore.' + side + '_'+prop][temp] = True
                self['throat.' + side + '_'+prop][tempt] = True
        if flowdir != 'all':
            temp2 =  flowdir
        else:
            temp2 = 'x'
        temp1 = ''
        if not replace:
            temp2 += '_0p' + str(ratio).split('.')[1]
            temp1 += '_0p' + str(ratio).split('.')[1]
        for elem in elements:
            self[elem + 'inlet'+temp1] = self[elem + 'inlet_'+temp2]
            self[elem + 'outlet'+temp1] = self[elem + 'outlet_'+temp2]

        #%% Save the results
        if replace:
            geo['pore.volume'][pinbc_id[offmask]] *= 0.5
            geo['pore.volume'][pinbc_id[offmask]] += geo['pore.volume'][pinbc_id[offmask]]*prat
            self.trim(pores=~pin)
            self._macro_Lx, self._macro_Ly, self._macro_Lz = Lnet2
            self._centpoint2 = cen
        else:
            tnew = self.find_neighbor_throats(pores=pin, mode='intersection')
            prop = 'all_net'+ endlabel
#            if suffix:
#                prop+= '_'+suffix
            self['pore.' + prop] = False
            self['throat.' + prop] = False
            self['pore.' + prop][pin] = True
#            self['pore.' + prop][pin] = True
            self['throat.' + prop][tnew] = True
            self._macro_Lx2, self._macro_Ly2, self._macro_Lz2 = Lnet2
            self._centpoint2 = cen
            self['prop.centpoint'+endlabel] = cen
            self['prop.Ldomain'+endlabel] = Lnet2

        self['pore.index'] = sp.arange(self.Np)
        self['throat.index'] = sp.arange(self.Nt)
        return

    def sphere_isect(self, center, radius=1e-6, pcoords=None):
        r'''
        Merge pores included in other pore volume.
        '''
        if pcoords is None:
            pcoords = self['pore.coords']

        distance = sp.power(sp.sum(sp.power(pcoords-center, 2), axis=1), 0.5)
        pmask = distance < radius
        return pmask


    def _dist_to_bc(self, p0, p1, psect, normal, eps=1e-12):
        r'''
        Find the intersectin point between p0-p1 line and plane identified
        by the normal line and a point in the plane (pplane).

         using vector eq: -----{ (x,y,z) = (x0,y0,z0) + t*m }-----
         with direction m = (x1,y1,z1) - (x0,y0,z0) = (m_x, m_y, m_z)
         (x0,y0,z0) and (x1,y1,z1) are the known coords of off- & in-bc pores
         and a scalar parametric t = (x-x0)/m_x = (y-y0)/m_y = (z-z0)/m_z
         if flowdir = 'x' --> t[left] = (x[left]-x0)/m_x[left]
                              t[right] = (x[right]-x0)/m_x[right]
         if flowdir = 'y' --> t[left] = (y[left]-y0)/m_y[left]
                              t[right] = (y[right]-y0)/m_y[right]
        '''
        nn = normal
        yy = psect - p0
        p0bc = p0 + yy.dot(normal)*normal
        yy = psect - p1
        p1bc = p1 + yy.dot(normal)*normal
        return p0bc, p1bc


    def _isect_line_plane(self, p0, p1, normal, pplane, eps=1e-12):
        r'''
        Find the intersectin point between p0-p1 line and plane identified
        by the normal line and a point in the plane (pplane).

         using vector eq: -----{ (x,y,z) = (x0,y0,z0) + t*m }-----
         with direction m = (x1,y1,z1) - (x0,y0,z0) = (m_x, m_y, m_z)
         (x0,y0,z0) and (x1,y1,z1) are the known coords of off- & in-bc pores
         and a scalar parametric t = (x-x0)/m_x = (y-y0)/m_y = (z-z0)/m_z
         if flowdir = 'x' --> t[left] = (x[left]-x0)/m_x[left]
                              t[right] = (x[right]-x0)/m_x[right]
         if flowdir = 'y' --> t[left] = (y[left]-y0)/m_y[left]
                              t[right] = (y[right]-y0)/m_y[right]
        '''
        linedir = p1 - p0
        ndotu = linedir.dot(normal)     # normal . linedir
        imask = sp.absolute(ndotu) > eps    # intersection mask
        if sp.all(~imask):
            msg = 'No line intersects the plane: normal = {} and point {}'
            print(msg.format(str(normal), str(pplane)))
        ww = p0[imask] - pplane     # p0 - pointPlane
        tf = -ww.dot(normal) / ndotu[imask]
        bcmask = sp.logical_and(tf >= 0, tf <= 1) # p0 - psect - p1
        inds = sp.where(imask)[0][bcmask]
        psect = p0[inds] + tf[bcmask].reshape((bcmask.sum(),1)) * linedir[inds]
        return psect, inds

    def center_point(self, pmin=None, pmax=None):
        r'''
        Calculate the center point of the network
        '''
        if pmin is None:
            pco = self['pore.coords']
            pmin = pco.min(axis=0)
            pmax = pco.max(axis=0)
        cen = sp.mean([pmin, pmax], axis=0)
        return cen


