# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:16:33 2015

@author: islah

===============================================================================
RegularLattice(Cubic): Generate lattice-like networks
===============================================================================

"""
import scipy as sp
import scipy.sparse as sprs
from numpy.linalg import norm
from OpenPNM.Network import Cubic
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class RegularLattice(Cubic):
    r"""
    This class generates a cubic network of the specified size and shape.
    Alternatively, an arbitrary domain shape defined by a supplied template.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : tuple of ints
        The (i,j,k) size and shape of the network.

    connectivity : int
        The number of connections to neighboring pores.  Connections are made
        symmetrically to any combination of face, edge or corners neighbors.

        Options are:

        - 6: Faces only
        - 8: Corners only
        - 12: Edges Only
        - 14: Faces and Corners
        - 18: Faces and Edges
        - 20: Edges and Corners
        - 26: Faces, Edges and Corners

    template : array of booleans
        An (i,j,k) array with True where the Network should be defined and
        False elsewhere. This approach is useful for creating networks of non-
        cuboid shape like spheres or cylinders, but still with a cubic lattice
        topology.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[3,4,5])
    >>> pn.Np
    60

    It is also possible to create Networks with cubic connectivity but
    non-Cubic shape by provding an array with True values where the network
    should exist to the ``template`` argument. The example below produces a sphere:

    >>> img = sp.ones([11, 11, 11])
    >>> img[5, 5, 5] = 0
    >>> from scipy.ndimage import distance_transform_bf as dt
    >>> img = dt(img) < 5  # Create a sphere of True
    >>> pn = OpenPNM.Network.Cubic(template=img)
    >>> pn.Np
    485

    If random distributions of coordination number is desired, one option is
    to create a Cubic network with many connections and the trim some:

    >>> pn = OpenPNM.Network.Cubic(shape=[5, 5, 5], connectivity=26)
    >>> Nt_original = pn.Nt
    >>> mod = OpenPNM.Network.models.pore_topology.reduce_coordination
    >>> pn.add_model(propname='throat.to_drop', model=mod, z=10, mode='random')
    >>> pn.trim(throats=pn['throat.to_drop'])
    >>> pn.Nt < Nt_original
    True
    """
    def __init__(self, shape=None, template=None, spacing=[1, 1, 1],
                 connectivity=6, **kwargs):
        super().__init__(shape, template, spacing, connectivity, **kwargs)
        self._connectivity = connectivity
        self['pore.inlet'] = self['pore.left']
        self['pore.outlet'] = self['pore.right']
        self._Lx = (shape[0]-1)*self._spacing[0]
        self._Ly = (shape[1]-1)*self._spacing[1]
        self._Lz = (shape[2]-1)*self._spacing[2]


    def add_boundaries(self, labels=['left', 'right']):
        r"""
        This method uses ``clone_pores`` to clone the surface pores (labeled
        'left','right', etc), then shifts them to the periphery of the domain,
        and gives them the label 'right_face', 'left_face', etc.
        """
        x, y, z = self['pore.coords'].T
        Lcx, Lcy, Lcz = self._spacing

        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0, 0, 0]
        offset['back'] = [Lcx*self._shape[0], 0, 0]
        offset['right'] = [0, Lcy*self._shape[1], 0]
        offset['top'] = [0, 0, Lcz*self._shape[2]]

        scale = {}
        scale['front'] = scale['back'] = [0, 1, 1]
        scale['left'] = scale['right'] = [1, 0, 1]
        scale['bottom'] = scale['top'] = [1, 1, 0]

        for label in labels:
            ps = self.pores(label)
            self.clone_pores(pores=ps, apply_label=[label+'_boundary', 'boundary'])
            # Translate cloned pores
            ind = self.pores(label+'_boundary')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords
        self['pore.index'] = sp.arange(0, len(self['pore.coords']))
        self['pore.inlet'] = self['pore.left_boundary']
        self['pore.outlet'] = self['pore.right_boundary']

    def trim_boundary_throats(self, inlet='inlet', outlet='outlet'):
        pinlet = self.pores(inlet)
        poutlet = self.pores(outlet)
        tinlet = self.find_neighbor_throats(pinlet)
        toutlet = self.find_neighbor_throats(poutlet)
        T = self['throat.conns']
        tbinlet = T[tinlet]
        indi = sp.array([throat[0] in pinlet and throat[1] in pinlet for throat in tbinlet])
        tboutlet = T[toutlet]
        indo = sp.array([throat[0] in poutlet and throat[1] in poutlet for throat in tboutlet])
        tout = sp.append(tinlet[indi], toutlet[indo])
        self.trim(throats=tout)

    def trim_direction_throats(self, directions=['z']):
        self.throat_wrt_direction()
        tdirects = self['throat.direction'] # x,y, z = 1, 2, 3
        adict = {'x':1, 'y':2, 'z':3}
        trim = [adict[way] for way in directions]
        tout = sp.zeros_like(tdirects, dtype=bool)
        for index in trim:
            tout += tdirects==index
        self.trim(throats=tout)

    def trim_row_throats(self, direction='z', rowstep=2, rowstart=1):
        '''Trim all throats in perpendicular planes to the supplied directions.
        For example, directions=['z'] and rowskip=1 mean trimming all throats in
        x-y planes and skip 1 x-y plane in between.
        '''
        adict = {'x':0, 'y':1, 'z':2}
        coords = self['pore.coords']
        conn = self['throat.conns']
#        row = conn[:, 0]
        col = conn[:, 1]
        dirpos = coords[:, adict[direction]]
        trows = dirpos[col]     # coords of throat's pore
        rows = sp.unique(dirpos) # unique coords of throat's pore
        tout = sp.zeros(self.Nt, dtype=bool)
        for row in rows[rowstart::rowstep]:
            ind = sp.where(trows==row)
            tout[ind] = True
        self.trim(throats=tout)

    def _to_direction(self, delta):
        r'''
        delta: coordinate difference between 2 pores of throat
        NOTE: opplicable for z = 26 only
        '''
        delta = sp.round_(delta)
        direction = sp.arange(1, 14)
        conditie = sp.array([[1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,0,-1], [0,1,-1], [-1,1,0], [1,1,1], [1,1,-1], [-1,1,1], [-1,1,-1]])
        result = sp.zeros(sp.shape(delta)[0])
        for i, no in enumerate(direction):
            ind = sp.where(sp.all(delta == conditie[i], axis=1))
            result[ind] = no
            ind = sp.where(sp.all(delta == -conditie[i], axis=1))
            result[ind] = no
        return result

    def create_connection_matrix(self, data=None, dropzeros=True,
                                 sym=True, sprsfmt='coo'):
        logger.debug('create_connection_matrix: Start of method')
        Np = self.num_pores()
        Nt = self.num_throats()

        # Check if provided data is valid
        if data is None:
            data = sp.ones((Nt,))
        elif sp.shape(data)[0] != Nt:
            raise Exception('Received dataset of incorrect length')

        # Clear any zero-weighted connections
        if dropzeros:
            ind = data > 0
        else:
            ind = sp.ones_like(data, dtype=bool)

        # Get connectivity info from network
        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        col = conn[:, 1]
        data_direction = self.throat_wrt_direction(ind)
        data = data*data_direction
        data = data[ind]

        # Append row & col to each other, and data to itself
        if sym:
            row = sp.append(row, conn[:, 1])
            col = sp.append(col, conn[:, 0])
            data = sp.append(data, data)

        # Generate sparse adjacency matrix in 'coo' format
        temp = sprs.coo_matrix((data, (row, col)), (Np, Np))

        # Convert to requested format
        if sprsfmt == 'coo':
            pass  # temp is already in coo format
        if sprsfmt == 'csr':
            temp = temp.tocsr()
        if sprsfmt == 'lil':
            temp = temp.tolil()
        logger.debug('create_connection_matrix: End of method')
        return temp

    def throat_wrt_direction(self, ind=None):
        r'''
        Represent throat in direction number (1-13)

        Ref: [1] a. Raoof and S. M. Hassanizadeh, “A new formulation for
        pore-network modeling of two-phase flow,” Water Resour. Res.,
        vol. 48, no. 1, pp. 1–13, 2012.
        '''
        if ind is None:
            ind = sp.ones(self.Nt, dtype=bool)
        conn = self['throat.conns'][ind]
        row = conn[:, 0]
        col = conn[:, 1]
        delta = (self['pore.coords'][col]-self['pore.coords'][row])/self._spacing
        self['throat.direction'] = self._to_direction(delta)
        return self['throat.direction']

    def randomize_mdpn(self, threshold=None):
        r'''
        Randomize network (construct an MDPNetwork) by labeling throat with
        either open or blocked state. The result is strored in
        self['throat.open'] in boolean array.

        Ref: Amir Raoof
        '''
        if threshold is None:
            threshold = sp.random.rand(13)
        T_thresh = self._to_threshold(threshold)
        elimination = sp.random.rand(self.Nt)
        self['throat.open'] = T_thresh > elimination

    def _to_threshold(self, threshold):
        r'''
        Represent throats in random threshold number
        '''
        direction = sp.arange(1, 14)
        try:
            tdirection = self['throat.direction']
        except:
            tdirection = self.throat_wrt_direction()

        result = sp.zeros_like(tdirection)
        for i, way in enumerate(direction):
            ind = sp.where(tdirection==way)
            result[ind] = threshold[i]
        return result

    def mdpn_match_connectivity(self, z, method='pso', lb=0, ub=1, kwargs={}):
        r'''
        Match the network's connectiviy to the specified
        mean (z=constant) or distribution (z=array) of connection number.

        z:   either a constant or an array-like quantity
        '''
        if isinstance(z, int):
#            self.match_connectivity_mean(z)
            threshold = [z/13]*13
            threshold = sp.array(threshold)
            self.randomize_mdpn(threshold)
            self.trim_blocked_throats()
        else:
            self.match_connection_matrix(z, method, lb, ub, kwargs)


    def match_connection_matrix(self, measured_z, method='pso', lb=0, ub=1,
                                 kwargs={}):
#                                maxiter=500, minfunc=0.5, nswarm=50, ):
        r'''
        Match the network's connectivity distribution to given measured
        connectivity distribution through manipulating the connection matrix.

        Arguments:
        ----------
        method:     'pso' for particle swarm optimization or
                    'ga' for genetic algorithm
        '''
        direction = sp.arange(1, 14)
        binz = sp.arange(1, 28)
#        z_target, binz = sp.histogram(measured_z, binz)
#        self._measured_connectivity = z_target
        conm = self.create_connection_matrix()

        kwargsout={'binz': binz, 'z_target': measured_z,
                   'conm': conm, 'direction': direction}
        if method == 'pso':
            try:
                from pyswarm.pso import pso
            except:
                logger.error('pyswarm is not installed')
            f_obj = self._randomize_connection_matrix      # objective function
            lb = [lb]*13     # lower bound
            ub = [ub]*13   # upper bound
            f_cons = None     # constraint function
            xopt, fopt = pso(f_obj, lb, ub, f_ieqcons=f_cons, kwargs=kwargsout,
                             **kwargs)
        elif method == 'pso_deap':
            kwargs.update(kwargsout)
            xopt, fopt = self._pso_deap(lb, ub, **kwargs)
            print('The best individu is: ', xopt,' : ', fopt)
        elif method == 'ga':
            kwargs.update(kwargsout)
            xopt, fopt = self._ga(lb, ub, **kwargs)
            print('The best individu is: ', xopt,' : ', fopt)

        self.randomize_mdpn(xopt)
        self.trim_blocked_throats()
        if kwargs['debug']:
            self._xopt = xopt
            self._fopt = fopt

    def _pso_deap(self, lb, ub, binz, z_target, conm, direction, omega=0.5,
                  npop=100, ngen=100, vlb=-1, vub=1, phip=0.5, phig=0.5,
                  **kwargs):
        r'''
        Arguments
        ---------
        - vlb/vub:  velocity lower/upper bounds
        '''
        from deap import base, creator, tools

        def generate(size, pmin, pmax, smin, smax):
            '''Initializer'''
            part = creator.Particle(sp.random.uniform(pmin, pmax, size))
            part.speed = sp.random.uniform(smin, smax, size)
            part.smin = smin
            part.smax = smax
            return part

        def updateParticle(part, best, omega, phi1, phi2):
            '''Updater'''
            u1 = sp.random.uniform(0, phi1, len(part))
            u2 = sp.random.uniform(0, phi2, len(part))
            v_u1 = u1 * (part.best - part)
            v_u2 = u2 * (best - part)
            part.speed *= omega
            part.speed += v_u1 + v_u2
#            for i, speed in enumerate(part.speed):
#                if speed < part.smin:
#                    part.speed[i] = part.smin
#                elif speed > part.smax:
#                    part.speed[i] = part.smax
            part += part.speed
            part[part<lb] = lb
            part[part>ub] = ub

        def f_obj(individu, binz, z_target, conm, direction):
            return (self._randomize_connection_matrix(
                    individu, binz, z_target, conm, direction),)

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Particle", sp.ndarray, fitness=creator.FitnessMin,
                       speed=list, smin=None, smax=None, best=None)
        toolbox = base.Toolbox()
        toolbox.register("particle", generate, size=13, pmin=lb, pmax=ub,
                         smin=vlb, smax=vub)
        toolbox.register("population", tools.initRepeat, list,
                         toolbox.particle)
        toolbox.register("update", updateParticle,
                         omega=omega, phi1=phip, phi2=phig)
        toolbox.register("evaluate", f_obj, binz=binz, z_target=z_target,
                         conm=conm, direction=direction)

        pop = toolbox.population(n=npop)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", sp.mean)
        stats.register("std", sp.std)
        stats.register("min", sp.amin)
        stats.register("max", sp.amax)

        logbook = tools.Logbook()
        logbook.header = ["gen", "evals"] + stats.fields

        best = None
        for g in range(ngen):
            for part in pop:
                part.fitness.values = toolbox.evaluate(part)
                if part.best is None or part.best.fitness < part.fitness:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
                if best is None or best.fitness < part.fitness:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values
            for part in pop:
                toolbox.update(part, best)

            # Gather all the fitnesses in one list and print the stats
            logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
            print(logbook.stream)

        f = lambda x: f_obj(x, binz=binz, z_target=z_target, conm=conm,
                            direction=direction)
        return best, f(best)[0]

    def _ga(self, lb, ub, binz, z_target, conm, direction, indpb=0.05,
            tournsize=3, npop=300, nhof=1, cxpb=0.5, mutpb=0.2, ngen=50,
            debug=False, **kwargs):
        r'''
        '''
        from deap import base, algorithms, creator, tools
        import random


        def threshold(xmin, xmax):
            x = random.random()
            x = xmin + x*(xmax-xmin)
            return x

        def cxTwoPointCopy(ind1, ind2):
            """Execute a two points crossover with copy on the input individuals. The
            copy is required because the slicing in numpy returns a view of the data,
            which leads to a self overwritting in the swap operation. It prevents
            ::

                >>> import numpy
                >>> a = numpy.array((1,2,3,4))
                >>> b = numpy.array((5.6.7.8))
                >>> a[1:3], b[1:3] = b[1:3], a[1:3]
                >>> print(a)
                [1 6 7 4]
                >>> print(b)
                [5 6 7 8]
            """
            size = len(ind1)
            cxpoint1 = random.randint(1, size)
            cxpoint2 = random.randint(1, size - 1)
            if cxpoint2 >= cxpoint1:
                cxpoint2 += 1
            else: # Swap the two cx points
                cxpoint1, cxpoint2 = cxpoint2, cxpoint1

            ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
                = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

            return ind1, ind2
        # Start GA
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))     # weight=-1 for minimization
        creator.create("Individual", sp.ndarray, fitness=creator.FitnessMin)

        toolbox = base.Toolbox()
        # Attribute generator
        #toolbox.register("attr_bool", random.randint, 0, 1)
        toolbox.register("attr_bool", threshold, lb, ub)
        # Structure initializers
        toolbox.register("individual", tools.initRepeat, creator.Individual,
                         toolbox.attr_bool, 13)     # individual of 13 thresholds
        toolbox.register("population", tools.initRepeat, list,
                         toolbox.individual)  # population of individual

        def f_obj(individu, binz, z_target, conm, direction):
            return (self._randomize_connection_matrix(
                    individu, binz, z_target, conm, direction),)

        toolbox.register("evaluate", f_obj,
                         binz=binz, z_target=z_target,
                         conm=conm, direction=direction)
        toolbox.register("mate", cxTwoPointCopy)
        toolbox.register("mutate", tools.mutFlipBit, indpb=indpb)
        toolbox.register("select", tools.selTournament, tournsize=tournsize)

        # def main_short():
#        random.seed(64)
        pop = toolbox.population(n=npop)
        hof = tools.HallOfFame(nhof, similar=sp.array_equal)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", sp.mean)
        stats.register("std", sp.std)
        stats.register("min", sp.amin)
        stats.register("max", sp.amax)

        pop, log = algorithms.eaSimple(pop, toolbox, cxpb=cxpb, mutpb=mutpb,
                                       ngen=ngen, stats=stats, halloffame=hof,
                                       verbose=debug)
        f = lambda x: self._randomize_connection_matrix(x, binz=binz,
                                                        z_target=z_target,
                                                        conm=conm,
                                                        direction=direction)
        fx = sp.array([f(item) for item in pop])
        ind = sp.where(fx==fx.min())[0][0]
        try:
            xopt = pop[ind]
            fopt = fx[ind]
        except:
            xopt, fopt = [pop[0], 1e10]
        return xopt, fopt

    def _randomize_connection_matrix(self, threshold, binz, z_target, conm,
                                     direction, seed=10, **kwargs):
        r'''
        Arguments:
        - threshold:    the parameter to be optimized
        - (binz, z):    connectivity distribution
        - conm:         connection matrix
        - direction:    the direction numbers
        '''
        conm_thresh = sprs.triu(conm)
        data = conm_thresh.data
        for i, direct in enumerate(direction):
            ind = sp.where(data == direct)
            data[ind] = threshold[i]
        if seed is not None:
            sp.random.seed(seed=seed)
        elimination = sp.random.rand(conm_thresh.nnz)   # uniform distribution
        Tmask_open = data > elimination

        z = sp.bincount(conm_thresh.row[Tmask_open])  # connection of each pore
        w_target = sp.ones_like(z_target)/len(z_target)
        w_z = sp.ones_like(z)/len(z)

        pdf_net, __ = sp.histogram(z, binz, weights=w_z)
        pdf, __ = sp.histogram(z_target, binz, weights=w_target)
        res = sp.sum(sp.absolute(pdf-pdf_net))
#        res = norm(pdf-pdf_net, ord=2)**2
        return res

    def trim_blocked_throats(self):
        r'''
        Trim all blocked throats after MDPN construction
        '''
        self.trim(throats=self.throats()[~self['throat.open']])

    def get_connectivity_mean(self, z_min=0):
        r'''
        mean connectivity for all throats with z > z_min
        '''
        T_nums = self.num_neighbors(self.pores())
        return sp.average(T_nums[T_nums > z_min])

# The following methods are from Original: OpenPNM.Network.models.pore_topology
    def match_connectivity_mean(self, z, trim=True):
        r"""
        Reduce the coordination number to the specified z value

        Parameters
        ----------
        z : int
            The (mean) coordination number or number of throats connected a pore

        Returns
        -------
        A label array indicating which throats are open and blocked to achieve
        desired coordination, stored in self['throat.open']. If trim=True, the
        blocked throats will be trimmed.

        Notes
        -----
        Pores with only 1 throat will be ignored in all calculations since these
        are generally boundary pores.

        """
        T_trim = ~self['throat.all']
        T_nums = self.num_neighbors(self.pores())
        # Find protected throats
        T_keep = self.find_neighbor_throats(pores=(T_nums == 1))

        z_ave = sp.average(T_nums[T_nums > 1])
        f_trim = (z_ave - z)/z_ave
        T_trim = sp.rand(self.Nt) < f_trim
        T_trim = T_trim*(~self.tomask(throats=T_keep))
        self['throat.open'] = ~T_trim
        if trim:
            self.trim_blocked_throats()

    def adjust_spacing(self, new_spacing, **kwargs):
        r"""
        Adjust the the pore-to-pore lattice spacing on a cubic network

        Parameters
        ----------
        new_spacing : float
            The new lattice spacing to apply

        Notes
        -----
        At present this method only applies a uniform spacing in all directions.
        This is a limiation of OpenPNM Cubic Networks in general, and not of the
        method.
        """
        coords = self['pore.coords']
        try:
            spacing = self._spacing
            self['pore.coords'] = coords/spacing*new_spacing
            self._spacing = new_spacing
        except:
            print('adjust spacing cannot be implemented')

    def get_subscripts(self, shape=None, **kwargs):
        r"""
        Return the 3D subscripts (i,j,k) into the cubic network

        Parameters
        ----------
        shape : list
            The (i,j,k) shape of the network in number of pores in each direction

        """
        if shape is None:
            shape = self._shape

        if self.num_pores('internal') != sp.prod(shape):
            print('Supplied shape does not match Network size, cannot proceed')
        else:
            template = sp.atleast_3d(sp.empty(shape))
            a = sp.indices(sp.shape(template))
            i = a[0].flatten()
            j = a[1].flatten()
            k = a[2].flatten()
            ind = sp.vstack((i, j, k)).T
            vals = sp.ones((self.Np, 3))*sp.nan
            vals[self.pores('internal')] = ind
            return vals
