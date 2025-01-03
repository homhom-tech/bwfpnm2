# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 16:39:23 2015

@author: islah
"""
import pickle
import os
from pandas import DataFrame, ExcelWriter
import scipy as sp
import psutil as pu
from bwfpnm.Utilities.__IO__ import *


def _make_geoinput(pradius, pvolume, pshapefactor, pconnectivity, pclayvolume,
                   tradius, tvolume, tlength, tshapefactor, tclayvolume,
                   tporelengths, tlengthtotal):
    r'''Collect all arguments into one dictionary'''

    geoinput = {'pradius': pradius,
                'pvolume': pvolume,
                'pshapefactor': pshapefactor,
                'pconnectivity': pconnectivity,
                'pclayvolume': pclayvolume,
                'tradius': tradius,
                'tvolume': tvolume,
                'tlength': tlength,
                'tshapefactor': tshapefactor,
                'tclayvolume': tclayvolume,
                'tporelengths': tporelengths,
                'tlengthtotal': tlengthtotal}
#                'tporeindices': tporeindices}
    return geoinput


def load_data(file_data, bcpores=True, netgeoinput=True):
    r'''Load a pickle data xxx.p file  and extract it into pores, throats,
    bc_throats, and macro.

    Arguments:
    ---------
    file_data:
    netgeoinput:    if True --> return (netinput, geoinput)
                    if False --> return (pores, throats, bc_throats, net[b'macro'])

    Example:
    ----------
    file_data = 'OpenPNM-develop/OpenPNM/Bwf/test/smallBenth.p'
    pores, throats, bc_throats, macro = load_data(file_data, False)
    or
    netinput, geoinput = load_data(file_data)
    '''
    with open(file_data, 'rb') as f:
        net = pickle.load(f, encoding='bytes')

    try:
        pores = net['pore']
        throats = net['throat']
        bc_throats = net['bc_throat']
    except:
        pores = net[b'pore']
        throats = net[b'throat']
        bc_throats = net[b'bc_throat']

    if bcpores:
        pass
#        _set_bcpores(pores, throats, bc_throats)

    if 'conns' not in throats.keys():
        conns = throats.pop('pores') - 1
        order = sp.argsort(conns, axis=1)
        grid = sp.ogrid[0:conns.shape[0], 0:2]
        conns = conns[grid[0], order]
        porelengths = throats.pop('porelengths')
        porelengths = porelengths[grid[0], order]
    else:
        conns = throats['conns']
        porelengths = throats['porelengths']

    if netgeoinput:
        netinput = {'coords': pores.pop('coords'),
                    'conns': conns,
                    'isInlet': pores.pop('isInlet'),
                    'isOutlet': pores.pop('isOutlet')}

        geoinput = _make_geoinput(pores.pop('radius'),
                                  pores.pop('volume'),
                                  pores.pop('shapefactor'),
                                  pores.pop('connectivity'),
                                  pores.pop('clayvolume'),
                                  throats.pop('radius'),
                                  throats.pop('volume'),
                                  throats.pop('length'),
                                  throats.pop('shapefactor'),
                                  throats.pop('clayvolume'),
                                  porelengths,
                                  throats.pop('lengthtotal'))
#                                  throats.pop('poreindices'))

        return (netinput, geoinput, net['macro'])
    else:
        return (pores, throats, bc_throats, net['macro'])


#==============================================================================
# def _set_bcpores(pores, throats, bc_throats):
# #    bc_throats: 'numbers', 'pores', 'radius', 'shapefactor', 'lengthtotal',
# #                'index', 'porelengths', 'length', 'volume', 'clayvolume'
#     # define the inlet-outlet direction to define the bcpore coordinates
#     pos_inlet = pores['coords'][pores['isInlet']]
#     pos_outlet = pores['coords'][pores['isOutlet']]
#     dpos = sp.absolute(pos_outlet - pos_inlet)
#     direction = sp.where(dpos==dpos.max())[0]    # flow direction: 0, 1, or 2
#     adir = sp.zeros(3, dtype=bool)
#     adir[direction] = True
#
#     # create new boundary pores from bc_throat['index']
#     nbc = bc_throats['numbers']
#     pcoords = sp.zeros((nbc, 3))      # bcpore coordinates
#     pconn = sp.ones(nbc)                        # bcpore connectvty no.
#     ppores = [[]]*nbc                   # bcpore neighbor pores
#     pisIn = sp.zeros(nbc, dtype=bool)
#     pisOut = sp.zeros(nbc, dtype=bool)
#     pthroats = [[]]*nbc
#     pindex = pores['index'][-1] + bc_throats['index']   # bcpore index
#     pvol = sp.zeros(nbc)
#     prad = sp.zeros(nbc)                        # bcpore radius
#     pG = sp.zeros(nbc)
#     pvolclay = sp.zeros(nbc)
#
#     props_t = ['radius', 'shapefactor', 'lengthtotal', 'porelengths',
#                'length', 'volume', 'clayvolume']
#     props_p = ['coords', 'connectivity', 'pore_neighbor', 'isInlet',
#                'inOutlet', 'throat_neighbor', 'index', 'volume', 'radius',
#                'shapefactor', 'clayvolume']
#
#     conns = bc_throats['pores']
#     for i, item in enumerate(conns):
#         mask = item == 0
#         mask += item == -1      # mask == True is the bc pore
#         pin = item[~mask]
#         # pore properties
#         pcoords[i] = pores['coords'][pin-1] # modify the x-coord with the throat length
# #        pcoords[i][adir] =
# #        pconn =
#         # throat properties
#         item[mask] = pindex[i]  # change pindex in the throat.conns
#
#         prad[i] = bc_throats['radius']
#==============================================================================


def _read_node1dat(filename=None):
    if filename is None:
        raise Exception('A file input is required!')

    with open(filename, 'r') as infile:
        # read first line
        num_pores, Lx, Ly, Lz = [float(x) for x in infile.readline().split()]
        num_pores = int(num_pores)
        domain_size = [Lx, Ly, Lz]

        index = sp.zeros(num_pores, dtype=int)
        pore_coords = sp.zeros([num_pores, 3])
        pore_connectivity = sp.zeros(num_pores, dtype=int)
        isInlet = sp.zeros(num_pores, dtype=bool)
        isOutlet = sp.zeros(num_pores, dtype=bool)

        pore_neighbor = [0]*num_pores
        pore_throat_neighbor = [0]*num_pores
        for i, line in enumerate(infile):   # read rest of lines
            # id, x, y, z, Z, idpore1, ..., idporen, isinlet, isoutlet,
            # idthroat1, ..., idthroatn
            array = line.split()
            index[i] = int(array[0])
            pore_coords[i] = [float(x) for x in array[1:4]]
            pore_connectivity[i] = int(array[4])
            neighbor_index = 4 + int(array[4]) + 1
            pore_neighbor[i] = [int(x) for x in array[5:neighbor_index]]
            isInlet[i] = int(array[neighbor_index])
            isOutlet[i] = int(array[neighbor_index+1])
            pore_throat_neighbor[i] = [int(x) for x in array[neighbor_index+2:]]

    return (num_pores, domain_size, pore_coords, pore_connectivity,
            pore_neighbor, isInlet, isOutlet, pore_throat_neighbor, index)


def _read_node2dat(num_pores, filename=None):
    if filename is None:
        raise Exception('A file input is required!')

    with open(filename, 'r') as infile:
        pore_volume, pore_radius = sp.zeros(num_pores), sp.zeros(num_pores)
        pore_shapefactor = sp.zeros(num_pores)
        pore_clayvolume = sp.zeros(num_pores)
        for i, line in enumerate(infile):   # read rest of lines
            # id, volume, radius, shapefactor/clayvolume
            array = [float(x) for x in line.split()]
            pore_volume[i], pore_radius[i] = array[1:3]
            pore_shapefactor[i], pore_clayvolume[i] = array[3:5]

    return (pore_volume, pore_radius, pore_shapefactor, pore_clayvolume)


def _read_link1dat(filename=None):
    if filename is None:
        raise Exception('A file input is required!')

    bound_pores = [-1, 0]       # -1: connected to inlet,   0: outlet
    with open(filename, 'r') as infile:
        # read first line
        num_throats = int(infile.readline())

        index = sp.zeros(num_throats, dtype=int)
        conns = sp.zeros([num_throats, 2], dtype=int)
        throat_radius = sp.zeros(num_throats)
        throat_shapefactor = sp.zeros(num_throats)
        throat_lengthtotal = sp.zeros(num_throats)

        bc_index = sp.zeros(num_throats, dtype=int)
        bc_conns = sp.zeros([num_throats, 2], dtype=int)
        bc_throat_radius = sp.zeros(num_throats)
        bc_throat_shapefactor = sp.zeros(num_throats)
        bc_throat_lengthtotal = sp.zeros(num_throats)
        i, j = 0, 0
        for line in infile:   # read rest of lines
            # Throat Index / Index of first connected pore/
            # Index of second connected pore/ Radius/ Shapefactor/ Total length
            array = line.split()
            pore1, pore2 = [int(x) for x in array[1:3]]
            if (pore1 not in bound_pores) and (pore2 not in bound_pores):
                index[i] = int(array[0])
                conns[i] = [pore1, pore2]
                throat_radius[i] = float(array[3])
                throat_shapefactor[i] = float(array[4])
                throat_lengthtotal[i] = float(array[5])
                i += 1
            else:
                bc_index[j] = int(array[0])
                bc_conns[j] = [pore1, pore2]
                bc_throat_radius[j] = float(array[3])
                bc_throat_shapefactor[j] = float(array[4])
                if pore1 == -1 or pore2 == -1:
                    # need to give sign for pn.add_inoutlet() method
                    bc_throat_lengthtotal[j] = -float(array[5])
                else:
                    bc_throat_lengthtotal[j] = float(array[5])
                j += 1
        num_throats = i
        bc_num_throats = j

    return (num_throats, conns[:i], throat_radius[:i], throat_shapefactor[:i],
            throat_lengthtotal[:i], index[:i],
            bc_num_throats, bc_conns[:j], bc_throat_radius[:j],
            bc_throat_shapefactor[:j], bc_throat_lengthtotal[:j], bc_index[:j])


def _read_link2dat(num_throats, bc_num_throats, filename=None):
    ''' Length of the first pore / length of the second pore /
    length of the throat / Volume / Clay volume'''
    if filename is None:
        raise Exception('A file input is required!')

    bound_pores = [-1, 0]
    with open(filename, 'r') as infile:
        throat_porelengths = sp.zeros([num_throats, 2])
#        throat_poreindices = sp.zeros([num_throats, 2])
        throat_length = sp.zeros(num_throats)
        throat_volume = sp.zeros(num_throats)
        throat_clayvolume = sp.zeros(num_throats)

        bc_throat_porelengths = sp.zeros([bc_num_throats, 2])
        bc_throat_length = sp.zeros(bc_num_throats)
        bc_throat_volume = sp.zeros(bc_num_throats)
        bc_throat_clayvolume = sp.zeros(bc_num_throats)
        i, j = 0, 0
        for line in infile:   # read rest of lines
            array = line.split()
            pore1, pore2 = [int(x) for x in array[1:3]]
            ## Read non boundary throats only
            if (pore1 not in bound_pores) and (pore2 not in bound_pores):
                throat_porelengths[i] = [float(array[3]), float(array[4])]
#                throat_poreindices[i] = [pore1, pore2]
                throat_length[i] = float(array[5])
                throat_volume[i] = float(array[6])
                throat_clayvolume[i] = float(array[7])
                i += 1
            else:   # BC throats only
                bc_throat_porelengths[j] = [float(array[3]), float(array[4])]
                bc_throat_length[j] = float(array[5])
                bc_throat_volume[j] = float(array[6])
                bc_throat_clayvolume[j] = float(array[7])
                j += 1

    return (throat_porelengths, throat_length,
            throat_volume, throat_clayvolume,
            bc_throat_porelengths, bc_throat_length,
            bc_throat_volume, bc_throat_clayvolume)


def _read_alldat(filelink='', filenames=None):
    if filenames is None:
        filenames = ['Berea_node1.dat', 'Berea_node2.dat', 'Berea_link1.dat',
                     'Berea_link2.dat']

    os.chdir(filelink)
    node1 = _read_node1dat(filenames[0])
    node2 = _read_node2dat(node1[0], filenames[1])
    link1 = _read_link1dat(filenames[2])
    link2 = _read_link2dat(link1[0], link1[6], filenames[3])
#    Orders of data:
#    ---------------
#    node1: 0-num_pores, 1-domain_size, 2-pore_coords, 3-pore_connectivity,
#           4-pore_neighbor, 5-isInlet, 6-isOutlet, 7-throat_neighbor, 8-throat_index
#    node2: 0-pore_volume,      1-pore_radius,
#           2-pore_shapefactor, 3-pore_clayvolume
#    link1: 0-num_throats,  1-conns[:i], 2-throat_radius[:i],
#           3-throat_shapefactor[:i], 4-throat_lengthtotal[:i], 5-index[:i],
#           6-bc_num_throats, 7-bc_conns[:j], 8-bc_throat_radius[:j],
#           9-bc_throat_shapefactor[:j], 10-bc_throat_lengthtotal[:j],
#           11-bc_index[:j]

#    link2: 0-throat_porelengths,       1-throat_length,
#           2-throat_volume,            3-throat_clayvolume,
#           4-bc_throat_porelengths,    5-bc_throat_length,
#           6-bc_throat_volume,         7-bc_throat_clayvolume

    macro = {'size': node1[1]}

    pore = {'numbers': node1[0],
            'coords': node1[2],
            'connectivity': node1[3],
            'pore_neighbor': node1[4],
            'isInlet': node1[5],
            'isOutlet': node1[6],
            'throat_neighbor': node1[7],
            'index': node1[8],
            'volume': node2[0],
            'radius': node2[1],
            'shapefactor': node2[2],
            'clayvolume': node2[3]}

    throat = {'numbers': link1[0],
              'pores': link1[1],
              'radius': link1[2],
              'shapefactor': link1[3],
              'lengthtotal': link1[4],
              'index': link1[5],
              'porelengths': link2[0],
              'length': link2[1],
              'volume': link2[2],
              'clayvolume': link2[3]}

    bc_throat = {'numbers': link1[6],
                 'pores': link1[7],
                 'radius': link1[8],
                 'shapefactor': link1[9],
                 'lengthtotal': link1[10],
                 'index': link1[11],
                 'porelengths': link2[4],
                 'length': link2[5],
                 'volume': link2[6],
                 'clayvolume': link2[7]}

    return (macro, pore, throat, bc_throat)


def _save_pickle(filename='net.p', **kwargs):
    with open(filename, 'wb') as outfile:
        pickle.dump(kwargs, outfile)


def save_dict_to_txt(filename='pnm.txt', **kwargs):
    with open(filename, 'w+') as outfile:
        for key in kwargs.keys():
            outfile.write(str(key) + '\t {}\n'.format(kwargs[key]))


def save_str_to_txt(string, filename='text.txt'):
    with open(filename, 'w+') as outfile:
        outfile.write(string)


def save_net_to_csv(net, label=True, prefix='berea'):
    r'''
    Save properties of a network and its related geometry to two csv files
    for pore and throat.
    '''
    for item in ['pore', 'throat']:
        keys = net.props(item)
        if label:
            keys.extend(net.labels(item))
        try:
            keys.remove(item+'.all')
        except:
            pass
        try:
            for geo in net._geometries:
                keys.extend(geo.props(item))
        except:
            pass
        table = {}
        for key in keys:
            if key not in ['pore.coords', 'throat.conns', 'throat.porelengths']:
                table[key] = net[key]
            elif key == 'pore.coords':
                table['x'] = net[key][:, 0]
                table['y'] = net[key][:, 1]
                table['z'] = net[key][:, 2]
            elif key == 'throat.conns':
                table['p1'] = net[key][:, 0]
                table['p2'] = net[key][:, 1]
            elif key == 'throat.porelengths':
                table['Lp1'] = net[key][:, 0]
                table['Lp2'] = net[key][:, 1]
        filename = prefix + '_' + item +'.csv'
        data = DataFrame(table)
        data.to_csv(filename, sep="\t")


def save_table_to_txt(table, filename='text.txt'):
    with open(filename, 'w+') as outfile:
        for row in table:
            try:
                string = ''.join(str(column).rjust(10) for column in row)
                outfile.write(string + '\n')
            except:
                outfile.write(str(row)+'\n')

def make_data_p(folder='SmallNetworkModel_PBModeled_Bentheimer_ECore/',
                name='Bentheimer1_smallNetwork',
                outfile=None, files=None):

    if outfile is None:
        outfile = name+'.p'

    if files is None:
        files = [name+'_node1.dat', name+'_node2.dat',
                 name+'_link1.dat', name+'_link2.dat']
    macro, pore, throat, bc_throat = _read_alldat(filelink=folder,
                                                  filenames=files)
    _save_pickle(filename=outfile, macro=macro, pore=pore, throat=throat,
                 bc_throat=bc_throat)


def _flow_pandas(obj):
    d = {}
    n = sp.size(obj._Pc)
    d['pc [Pa]'] = obj._Pc
    d['log(pc)'] = obj._log_Pc
    d['saturation [-]'] = obj._saturation
    d['moisture content [kg/m3]'] = obj._moisturecontent

    try:
        if sp.size(obj._saturation_surf) == n:
            d['saturation surface [-]'] = obj._saturation_surf
        else:
            d['saturation surface [-]'] = sp.zeros(n)
    except:
        d['saturation surface [-]'] = sp.zeros(n)

    try:
        if sp.size(obj._saturation_vapour) == n:
            d['saturation vapour [-]'] = obj._saturation_vapour
        else:
            d['saturation vapour [-]'] = sp.zeros(n)
    except:
        d['saturation vapour [-]'] = sp.zeros(n)

    try:
        if sp.size(obj._permeability) == n:
            d['permeability effective [s]'] = obj._permeability
        else:
            d['permeability effective [s]'] = sp.zeros(n)
    except:
        d['permeability effective [s]'] = sp.zeros(n)
    try:
        d['permeability relative [-]'] = obj._rel_permeability
    except:
        d['permeability relative [-]'] = sp.zeros(n)
    try:
        d['permeability absolute [m2]'] = obj._abs_m2_permeability
    except:
        d['permeability absolute [m2]'] = sp.zeros(n)

    try:
        d['permeability absolute [mD]'] = obj._abs_mD_permeability
    except:
        d['permeability absolute [mD]'] = sp.zeros(n)

    try:
        d['conductance effective [sm]'] = obj._conductance
    except:
        d['conductance effective [sm]'] = sp.zeros(n)

    df = DataFrame(data=d)
    return df


def _data_pandas(data):
    r'''
    Data argument must be a dictionary with an array/list value per each key.
    The arrays all must have the same size.
    '''
    d = {}
    for key, val in data.items():
        d[key] = val

    df = DataFrame(data=d)
    return df


def save_flow_csv(obj, filename=''):
    r'''
    Save all important data of flow algorithm (moisture content, permeability and related variables) into a csv file.

    Arguments:
    ----------
    obj:    an instance of flow object
    filename:   a string of file name. a file type .csv will be added automatically if not pprovided

    '''
    if filename == '':
        filename = obj.name + '.csv'
    else:
        filename = filename.split('.')[0] + '.csv'

    df = _flow_pandas(obj)
    df.to_csv(filename, sep="\t")


def save_data_csv(data, filename=''):
    r'''
    data must be a dictionary with name and its value
    '''
    if filename == '':
        filename =  'data.csv'
    else:
        namesplit = filename.split('.')
        if namesplit[-1] != 'csv':
            filename += '.csv'

    df = _data_pandas(data)
    df.to_csv(filename, sep="\t")


def read_flow_csv(filename=''):
    r'''
    Read a csv saved from an instance of flow algorithm (as a result of permeability calculation). This function is the inverse of 'save_flow_csv' function.

    Return (header, pc, lpc, w, k)
    '''
    pc, lpc, w, k = [], [], [], []

    with open(filename, 'r') as infile:
        # read first line
        header = [x for x in infile.readline().split('\t')]

        for line in infile:   # read rest of lines
            # no, log(pc), moisture content [kg/m3],
            # pc [Pa], Kabs [m2], Kabs [mD], Keff [s], Krel, S, Ssurf, Svap
            array = [float(x) for x in line.split()]
            pc.append(array[3])
            lpc.append(array[1])
            w.append(array[2])
            k.append(array[6])

    return (header, pc, lpc, w, k)


def save_sparse_csr(filename, csr_matrix):
    sp.savez(filename, data=csr_matrix.data , indices=csr_matrix.indices,
             indptr=csr_matrix.indptr, shape=csr_matrix.shape )


def load_sparse_csr(filename):
    loader = sp.load(filename)
    csr = sp.sparse.csr_matrix
    return csr((loader['data'], loader['indices'], loader['indptr']),
               shape=loader['shape'])


def save_cpu_log():
    pid = pu.os.getpid()
    print("%CPU\t%MEM")
    try:
        while True:
            x = _get_cpumem(pid)
            if not x:
                print("no such process")
                exit(1)
            print("%.2f\t%.2f" % x)
            time.sleep(1)
    except KeyboardInterrupt:
        print
        exit(0)


def _get_cpumem(pid=None):
    if pid is None:
        pid = pu.subprocess.os.getpid()

    d = [i for i in pu.subprocess.getoutput("ps aux").split("\n")
        if i.split()[1] == str(pid)]
    return (float(d[0].split()[2]), float(d[0].split()[3])) if d else None


def load_properties(file_data=None):
    r'''Load the fluid property data file
    '''
    if file_data is None:
        file_data = '/bwfpnm/Phase/bwfpnm_parameters.dat'

    props = {}
    with open(file_data, 'r') as infile:
        for row in infile:
            row = row.split()
            if not row:
                continue
            if row[0] == 'phase':
                props[row[1]] = {}
                phase = row[1]
            elif row[1].strip() != 'None':
                props[phase][row[0]] = float(row[1])
            else:
                props[phase][row[0]] = None
    return props


def read_raw(self, filename='')  :
    '''Read a raw file with scipy and return it as a matrix
    '''
    A = sp.fromfile(filename, dtype='int16', sep="")
#    A = A.reshape([1024, 1024])
#    plt.imshow(A)
    return A



#def save_flow_xls(obj, filename=''):
#    if filename == '':
#        filename = obj.name + '.xls'
#    else:
#        filename = filename.split('.')[0] + '.xls'
#
#    df_out = _output_pandas(obj)
#    writer = ExcelWriter(filename)
#
#    df_out.to_excel(writer, sheet_name='flow', engine='xlwt')
#    writer.save()
