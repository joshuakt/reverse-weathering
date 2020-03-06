Version 1.0

This set of python scripts runs our sediment diagenesis model, as described in Krissansen-Totton and Catling (2020) "A coupled carbon-silicon cycle model over Earth history: Reverse weathering as a possible explanation of a warm mid-Proterozoic climate", Earth and Planetary Science Letters. The outputs from this sediment diagenesis model are then fit with polynomials, and the resulting parameterizations are used in our coupled carbon-silica cycle model (github.com/joshuakt/reverse-weathering/carbon-cycle-with-RW)

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton and Catling (2020). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (jkt@ucsc.edu)

REQUIREMENTS: Python, including numpy, pylab, and scipy modules.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run_sediment_diagenesis_model.py 
This script runs the sediment diagenesis model for a range of ocean boundary conditions i.e. different dissolved silica concentrations. This can be done for all the scenarios in the main text (Fig. 3 - 7). Additionally, the sediment diagenesis model can be run either with or without biogenic silica precipitation, representing cases in the Phanerozoic and Precambrian, respectively. Outputs are save as numpy arrays.

%% Calculate_fluxes_func.py
Sediment diagenesis model for case with biogenic silica precipitation (Phanerozoic)

%% Calculate_fluxes_func_Aonly.py
Sediment diagenesis model for case withot biogenic silica precipitation (Precambrian)

%% parameter_fit_to_diagenesis_model.py
This script loads the outputs from the sediment diagenesis models, and fits them with simple functions (usually polynomials) for use in the coupled carbon-silica cycle model. Specifically, these are used in the script Si_cycle_functions.py. Polynomial functions are chosen heuristically.

%% mcmc_parameter_fit.py
This script uses a MCMC approach to find sediment diagenesis parameters that fit modern conditions. Code needs to be run separately for proximal, distal, and pelagic conditions. Choice of T&D modern RW flux and Rahman et al. modern RW flux.

-------------------------
Contact e-mail: jkt@ucsc.edu
