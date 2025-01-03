BwfPNM
=======

.. contents::


Installation
------------
Summary:
--------
'''
conda create -n envname python matplotlib pandas pyamg spyder scikit-umfpack
pip install pyswarm OpenPNM==1.3
pip install git+https://islah@bitbucket.org/islah/bwfpnm.git
'''

Detail:
----------
1. Download and install Anaconda, please go to its webpage.

2. Update and configure an anaconda environment. In a terminal do the following:
	a. update conda: conda update anaconda
	b. update spyder: conda update spyder
	c. create an environment: conda create -n [envname]
	d. activate the env: source (CONDA) activate [envname]

2. [Optional] Install scikit-umfpack
	conda install -c conda-forge scikit-umfpack, or
	conda install -c zeehio scikit-umfpack, or
	pip install scikit-umfpack

2. Install pyswarm and pyamg
	First, try the following:
		conda install pyamg
		conda install -c auto pyswarm (LUKTE NIET)
	Use pip if that doesn't work
		pip install pyswarm pyamg (DIT WEL MET ENKEL PYSWARM)

2. Install numpy, scipy, matplotlib, pandas, spyder
	conda install scipy numpy matplotlib pandas spyder==3.2	(SPYDER 3.2 NIET GELUKT)
3. Install OpenPNM 1.3

    	pip install OpenPNM==1.3
	
	Note: you may be asked to install matplotlib, dateutil.
	If error occured due to dateutil-related problem, try:
	pip install matplotlib

4. Replace the installed OpenPNM folder with the supplied OpenPNM folder as it has been modified. (GEEN IDEE) De folder zoeken en ergens opslaan -> path intypen in explorer
	- linux: /home/anaconda/env
	- windows: /users/[user name]/anaconda/env

5. Install bwfpnm

	pip install -e "[path to bwfpnm/]"

5.a. install a specific commit:
	pip install git+https://islah@bitbucket.org/islah/bwfpnm.git@2418c9db5ca2014b13e4e6bd7595fd6adb507e18 (NIET GELUKT)



The following modules may also be installed in the environment.

Required modules
----------------
- anaconda > 4.1.1
- python > 3.5.2
- scikit-umfpack
- scipy > 0.18.0
- numpy > 1.11.1
- pandas > 0.18
- matplotlib
  - installing dateutil may result in an error. Do the following instead: 
'''
conda install python-dateutil
'''
- spyder > 2.3.9
- openpnm 1.3
  pip install OpenPNM==1.3
- bwfpnm

Optional
---------------
- scikit-umfpack
	conda install -c conda-forge scikit-umfpack, or
	conda install -c zeehio scikit-umfpack, or
	pip install scikit-umfpack
- pyamg 3.0.2
- pyswarm 0.6


Simulation
------------
1. Convert the four '.dat' files into a '.p' file (convert_data.py)

2. Run the desired invasion algorithm (invasion.py)

3. Calculate the hygric properties (permeability.py)

4. Others:
	a. Network properties (net_properties.py)
	b. Create a paraview file: ctrl.export(network, name)
