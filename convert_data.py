#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:07:36 2018

@author: islah

Two multiscale networks from Jiang: mulnet & mulnet_big. They have the same
format as that of Statoil because they do not contain hydraulic radius data.

They are two ways to convert them:
    1. Trimming the supplied boundary throats (and pores).
        --> bpnm.Utilities.IO.make_data_p(folder, name, outfile)

    2. Modifying the supplied boundary throats such that they are located at
        the inlet/outlet plane while preserving the geometrical properties.
        --> pn = bpnm.Utilities.IOnew.Statoil.load(path, prefix, suffix)

In the future, the two-mentioned functions should be merged.
"""

import bwfpnm as bpnm
import os  # This module provides a portable way of using operating system dependent functionality


# loc = os.path.dirname(os.path.dirname(loc))
for i in range(1):
    loc = os.getcwd()
    name = 'sintered_glass'  # name of the pore network to be converted
    # name = 'CB_7mm_1.2um_1_1800.raw'
    bc = True
    csv_net = False
    # %% Convert data *.dat files to a *.p file
    loc = os.path.join(loc, name)
    if not bc:
        # %% Convert data *.dat files to a *.p file by trimming the boundary throats
        subdir = loc
        bpnm.Utilities.IO.make_data_p(folder=subdir, name=name, outfile=loc + '.p')
    else:
        # %% Convert *.DAT files with additional boundary pores and throats
        # --------------------------------------------------------------------------
        path = loc
        pn = bpnm.Utilities.IOnew.Statoil.load(path=path, prefix=name, suffix='BC')
        if csv_net:
            bpnm.Utilities.IO.save_net_to_csv(pn, label=True, prefix=name.lower() + '_BC')

