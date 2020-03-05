Version 1.0

This set of python scripts runs our geological carbon cycle model, now updated to include reverse weathering, for the last 4.0 Ga and plots selected outputs alongside proxy data from the literature. Python code for the sediment diagenesis model is also provided. Both models are described in Krissansen-Totton and Catling (2020) "A coupled carbon-silicon cycle model over Earth history: Reverse weathering as a possible explanation of a warm mid-Proterozoic climate", Earth and Planetary Science Letters. The carbon cycle code is adapted from Krissansen-Totton et al. (2018; PNAS), the original code for which is also available at github.com/joshuakt/early-earth-carbon-cycle.

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton and Catling (2020). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (jkt@ucsc.edu)

REQUIREMENTS: Python, including numpy, pylab, and scipy modules.

HOW TO RUN CODE:
(1) Put all the python scripts in the same directory, and ensure python is working in this directory.
(2) Open time_series_parallel.py and check desired parameter ranges, options, and number of iterations (default parameters reproduce Fig. 3 in main text)
(3) Run time_series_parallel.py. Code will output model confidence interval alongside proxy data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF CODE STRUCTURE:

%% time_series_parallel.py
This script provides the shell to repeatedly call the forward model and plot the output distributions. The "Options" section of the code contains various parameters that can be modified by the user:

- it_num: number of forward model calls used to build output distributions. The larger this is, the longer the code will take to run.
- Parallelize: determines whether multiple threads should be used (faster) or whether all iterations should be computer using the same thread (slower).
- Carbon_chem: determines whether carbon equilibrium chemistry constants should be temperature dependent (slower but more accurate) or fixed (faster but less accurate).
- methane: controls whether methane should be added to the Proterozoic and Archean atmospheres

Ranges for uncertain parameters can also be modified by the user. Modifying anything else in 'time_series_parallel.py' may produce errors in the code. Once parameter ranges have been defined, the script calls the forward model ('Forward_Model') once and uses the dimensions of the outputs to define an output array to store all future outputs ('all_output'). The forward model is then called many times (equal to 'it_num') and parameter ranges are randomly sampled for each forward model call. Forward model calls resulting in errors or non-physical outputs are discarded (the code may still print error messages from these discarded outputs). The remainder of the script calculates 95% confidence intervals for model outputs based on the distribution of outputs, and plots these alongside proxy data and other literature estimates (see below).

%% time_series_functions_reformate.py
This script contains the following functions which, taken together, define and solve the forward model:

% Forward_Model - Given parameter inputs, Forward_Model calculates the initial (modern) conditions for the carbon cycle e.g. equilibrium ocean chemistry and modern fluxes. Proportionality constants for carbon cycle functions are also calculated from initial conditions. Next, the ODE solver is called, and the system of equations describing the carbon cycle are solved. The ODE solver only returns DIC and ALK as a function of time for both the ocean and the pore space. These outputs are fed back into carbon_cycle to obtain the time evolution for carbon cycle fluxes and ocean chemistry. Selected outputs are returned to time_series_parallel.py. 

% system_of_equations - Contains the ODEs that describe the time evolution of the carbon and silica cycles. The function takes the current state of the carbon and silica cycle (calculated using the carbon_cycle function), and returns the time derivatives of DIC and ALK in the ocean and the pore space, as well as dissolved silica concentrations. This function is fed into the ODE solver to compute the time evolution of DIC, ALK for the ocean and pore space, and [Si] for the ocean.

% lf - Calculates continental land fraction as a function of time.

% carbon_cycle - This function computes equilibrium chemistry and carbon cycle fluxes from DIC and ALK for both the ocean and the pore space. First, carbon_cycle calculates equilibrium chemistry for the ocean and the pore space, then heatflow and outgassing parameters. Next, it calculates surface and pore space temperatures from pCO2 using our climate model and surface-deep ocean relationship. Note that the option 'Carbon_chem' in time_series_parallel.py determines whether this procedure is done once (fast, less accurate), or several times such that the correct temprature is converged upon for calculating carbon chemistry equilibrium constants (slow, more accuate). Finally, carbon cycle fluxes are calculated from this information. Both carbon cycle fluxes and equilibrium chemistry are returned to system_of_equations or Forward_Model.

%% thermodynamic_variables.py
Contains the function Sol_prod that calculates carbonate solubility product as a function of temperature from Pilson, M. E. (1998) "An Introduction to the Chemistry of the Sea", Prentice-Hall, Inc. This function reproduces table G.1 in Pilson for salinity=35. thermodynamic_variables.py also contains the temperature-dependent functions for carbon equilibrium constants, equil_cont.

%% clim_model.py
This contains the the parameterized climate models used in this study. clim_fun_CO2_only is the polynomial fit to the CO2-only climate model, clim_fun_lowCH4 is the polynomial fit to the 100 ppm methane model, and clim_fun_highCH4 is the polynomial fit to the 1% methane model.

%% Si_cycle_functions.py
Contains functions for abiotic silica precipitation, biotic silica precipitation, and reverse weathering. These parameterizations are derived from the sediment diagenesis model. Different cases in the main text (Fig. 3, 4, 5, 6, 7) require different choices of parameterizations for reverse weathering and biotic silica functions.

% RW_flux - Parameterization of revererse weathering silica sink

% Si_Bio_flux - Parameterizations of biogenic silica sink

% Si_aBio_flux - Parameterizations of abiotic silica sink

%% Other files included in the code:
- Halevy_Bachan_high.txt and Halevy_Bachan_low.txt are ocean pH estimates from Halevy and Bachan than are plotted along our model outputs for compairson
- options_array.npy is an array created by the code for storing user choices about parallelization, methane, iteration number, and carbon chemistry precision.
- original_XXX.npy - nominal model outputs from Krissansen-Totton et al. (2018) to compare to no reverse weathering case
- Conley_data.txt - Estimates of ocean dissolved silica through time from Conley et al.
END EXPLANATION OF CODE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

-------------------------
Contact e-mail: jkt@ucsc.edu
