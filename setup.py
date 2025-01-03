#!/usr/bin/env python

import os
import sys
sys.path.append(os.getcwd())

import scipy as sp

if sp.__version__ < '0.14.0':
	raise Exception('bwfpnm requires SciPy version 0.14.0 or greater')

try:
    import OpenPNM as pnm
except:
    raise Exception('bwfpnm requires OpenPNM version 1.1-beta or greater')

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='bwfpnm',
    packages=['bwfpnm',
              'bwfpnm.Base',
              'bwfpnm.Network',
              'bwfpnm.Geometry',
              'bwfpnm.Geometry.models',
              'bwfpnm.Phases',
              'bwfpnm.Phases.models',
              'bwfpnm.Physics',
              'bwfpnm.Physics.models',
              'bwfpnm.Utilities',
              'bwfpnm.Algorithms',
              'bwfpnm.Postprocessing'],
    version='1.0-beta',
    description="An object-oriented code for simulating pore network model specifically for porous building materials using OpenPNM framework.",
    author='BwfPNM Team: Islah',
    author_email='m.islah22@gmail.com',
    download_url='https://github.com/islah/BwfPNM/',
    url='https://github.com/islah/BwfPNM',
	install_requires = ['scipy>=0.14.0', 'OpenPNM>=1.1-beta'],
)
