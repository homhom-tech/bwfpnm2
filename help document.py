# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 09:30:28 2019

@author: u0131531
"""

import bwfpnm as bpnm
import scipy as sp
from numpy.linalg import cond

ctrl = bpnm.Base.Controller()
ctrl.loglevel = 40

file_data = 'OpenPNM-develop/OpenPNM/Bwf/test/smallBenth.p'
pores, throats, bc_throats, macro = bpnm.Utilities.IO.load_data(file_data, False)