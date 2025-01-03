# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 18:58:08 2015

@author: islah
"""
from OpenPNM.Geometry import models as gm
#from OpenPNM.Bwf.Geometry import models as bgm
from bwfpnm.Geometry import GenericGeometry
#from .models import __poredata__ as poredata


class Throat(GenericGeometry):
    r"""
    Throat geometry for throat model.
    Pore's diameter, volume, area, and length are set to zero.

    Parameters
    ----------
    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.
    tdiameter, tlength, tarea:  throat's diameter, length, and area

    """

    def __init__(self, tdiameter, tlength, tarea=None, **kwargs):
        r"""
        Initialize
        """
        super(Throat,self).__init__(**kwargs)
        self._generate(tdiameter, tlength, tarea)

    def _generate(self, tdiameter, tlength, tarea):
        self['pore.diameter'] = 0.0
        self['pore.volume'] = 0.0
        self['pore.area'] = 0.0
        self['pore.length'] = 0.0

#        if n_data is None:
#            n_data = self.Nt
#        rand_data, poredistribution = poredata.paperbased_radii_2d(
#            material, case, n_data, fig_plot=False)

        self['throat.diameter'] = tdiameter
        self['throat.length'] = tlength
        try:
            self.models.add(propname='throat.volume',
                            model=gm.throat_volume.cylinder)
            self.models.add(propname='throat.area',
                            model=gm.throat_area.cylinder)
        except:
            pass
        if tarea is not None:
            self['throat.area'] = tarea


if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
