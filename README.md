# Overview

*bwfpnm* is a software to determine the moisture storage and transport properties of porous building materials from their pore-network models. This software is a customization of the open-source [OpenPNM 1.3](https://github.com/PMEAL/OpenPNM/releases/tag/V1.3).

## Installation

### Summary
For fast installation, please run the following:
```
conda create -n [sim] python==3.6 matplotlib pandas pyamg spyder scikit-umfpack
source activate [sim]
pip install pyswarm OpenPNM==1.3
pip install git+https://islah@bitbucket.org/islah/bwfpnm.git
```

If you encounter an error, please do the detail steps instead.


### Detail
1. Download and install Anaconda, please go to its webpage.

2. Update and configure an anaconda environment. In a (conda) terminal do the following:
    
    * update conda: `conda update anaconda`
    * update spyder: `conda update spyder`
    * create an environment: `conda create -n [envname]`
    * activate the env: `source activate [envname]`

3. (Optional) Install scikit-umfpack:

	`conda install -c conda-forge scikit-umfpack`, or

	`conda install -c zeehio scikit-umfpack`, or

	`pip install scikit-umfpack`

4. (Optional) Install pyswarm and pyamg.

    First, try the following:

    `conda install pyamg`
    `conda install -c auto pyswarm`

    Use pip if the above doesn't work

    `pip install pyswarm pyamg`

5. Install numpy, scipy, matplotlib, pandas, spyder

    `conda install scipy numpy matplotlib pandas spyder`

6. Install OpenPNM 1.3

    `pip install OpenPNM==1.3`
	
    Note: you may be asked to install matplotlib, dateutil. If error occured due to dateutil-related problem, try:
    `pip install matplotlib`

7. Replace the installed OpenPNM folder with the supplied customized OpenPNM folder. You may find the location at:
    - linux: /home/anaconda/env/lib/python3.6/site-package
    - windows: /users/[user name]/anaconda/env/...

8. Install bwfpnm

    `pip install -e "[path to bwfpnm/]" `

    a. To install a specific commit instead:

        `pip install git+https://islah@bitbucket.org/islah/bwfpnm.git@2418c9db5ca2014b13e4e6bd7595fd6adb507e18`


## Example Usage

The files related to this section can be found in the folder `working_example`. Other examples may be outdated.

In general, there are three steps to do a simulation:

1. Convert the four `.dat` files of the Statoil format to a single file of native python format `.p`.

    * Please refer to the sample file `convert_data.py` for further detail

2. Run an invasion algorithm to determine the moisture distribution at each condition/capillary pressure.
    * Please refer to the sample file `invasion.py` for further detail. Here are the steps done in the file:
        1. Load the network data in the '.p' file.
        2. Create all necessary objects: network, geometry, phases, and their physics. And define all parameters required for the simulation.
        3. Run the invasion algorithm to create the invasion object(s).
        4. Save the invasion object(s) to a '.pnm' file.

3. Run the permeability algorithm to determine the moisture retention and permeability curves.
    * Please refer to the sample file `permeability.py` for further detail. Here are the steps done in the file:
        1. Load the invasion object(s) from the '.pnm' file.
        2. Define all parameters required for the simulation.
        3. Run the permeability algorithm to create the permeability objects.
        4. Save the results into .csv files. Note that the permeability object cannot be saved.


The above steps are based on a module containing all routine steps for a moisture flow simulation. This module is called `bwfr` in the sample files. You can always make your own version if desired.

There is the file `bwfpnm_parameters.dat` containing all parameters that are probably not mentioned in the script files. Have a look and alter it as desired.

Prior to the above steps, the environment must be first activated:
`source activate [envname]`
in Windows, do this instead:
`activate [envname]`
then open the editor (spyder):
`spyder`

