import numpy as np
import pylab

## This file contains the functions for reverse weathering and biogenic silica fluxes, as determined by polynomial fits to the outputs of the sediment diagenesis model.
## The different fits represent different parameter choices, as described in the manuscript. Functions can be commented out to try out different parameter choices 
## and reproduce the different figures in the main text.

def RW_flux(Si_abun,H_ocean,t0,initRW):
    Si_abun = Si_abun*1000 # convert mMol/kg from mol/kg

##############################################
#############################################

### RW parameterization for Fig. 3 in main text
    f1 = 60.38939292*Si_abun**5 -126.21506006 *Si_abun**4 + 69.16150651*Si_abun**3  + 3.33014549* Si_abun**2  + 1.38272369*Si_abun + 0.65007329
    if Si_abun< 0.186472:
        f2 = 0.0
    else:
        f2 = np.max([0.0,2.03098274*Si_abun**6 -16.96872155*Si_abun**5 + 54.66173384*Si_abun**4  -84.21016421* Si_abun**3  + 61.24035109*Si_abun**2 - 9.86063672*Si_abun + 0.18330492])

### RW parameterization for Fig. 4 in main text
#    f1 =  Si_abun**4*46.38053287+Si_abun**3*-91.02008426  +Si_abun**2 *51.72188425 + Si_abun *  0.33731983  + 4.041267
#    if Si_abun < 0.0595411:
#        f2 = 0.0
#    else:
#        f2 = 0.06543753*Si_abun**5  -0.8972411 *Si_abun**4 + 3.92389598*Si_abun**3  -7.39240127* Si_abun**2  + 10.28661551*Si_abun -0.58318515

#### RW parameterization for Fig. 5 in main text
#    f1 =   1.08847555 *np.tanh(Si_abun*344.37396583-223.93762487) + 0.83950234 * Si_abun + 2.50457139
#    f2 = np.max([0,1.17084241 *np.tanh(Si_abun* 163.39926589 - 111.42260266) + 0.74820693 * Si_abun + 1.01128036])

### RW parameterization for Fig. 6 in main text
#    f1 = -34.68555321*Si_abun**3 +  77.32055905*Si_abun**2 -4.54075261* Si_abun  +  1.75363474
#    f2 = np.max([0.0,  -6.51803015* Si_abun**3  + 26.1737101*Si_abun**2 +15.37160264*Si_abun -3.02992175])
 

### RW parameterization for Fig. 7 in main text
#    f1 = 3.6166122*Si_abun**3 + 0.30710427*Si_abun**2 + 4.83648973*Si_abun + 1.15605322 
#    f2 = np.max([0.0, -3.01022775 * Si_abun**3  + 14.64895684*Si_abun**2 -8.39560695 *Si_abun + 0.6553986]) 
#    if Si_abun < 0.5:
#        f2 = 0.0 
    
    ## Gradual transition between Phanerozoic and Precambrian parameterizations:
    if t0 <= 0.541e9:
        flux = f1*1e12
    elif (t0>0.541e9)and(t0<1.041e9):
        flux = 1e12* (f1*(1.041e9-t0)/5e8  + f2*(t0-0.541e9)/5e8)
    else:
        flux = f2*1e12 
    
    if Si_abun<0:
        flux =0.0

    return flux


def Si_Bio_flux(Si_abun,t0): 
    Si_abun = Si_abun*1000 # convert mMol/kg from mol/kg

##############################################
#############################################

### Biogenic silica parameterization for Fig. 3 in main text
    f1 = 1725.4786102*Si_abun**6  -3541.63689444*Si_abun**5 + 2822.12738589*Si_abun**4  - 1039.66772491* Si_abun**3  + 188.47492215*Si_abun**2 -4.04821308*Si_abun + 8.00302464
    f2 = 0.0*Si_abun

## Biogenic silica parameterization for Fig. 4 in main text
#    f1 = 4.22121400e+03 *Si_abun**6  -7.99808174e+03*Si_abun**5 +  5.87221586e+03*Si_abun**4  -1.97706858e+03* Si_abun**3  + 3.19644569e+02 *Si_abun**2 -7.55719889e+00*Si_abun + 6.92737494e+00
#    f2 = 0.0*Si_abun

#### Biogenic silica parameterization for Fig. 5 in main text
#    f1 =  6.13182221e+03*Si_abun**6  -1.14883069e+04*Si_abun**5 +  8.33995255e+03*Si_abun**4   -2.76295885e+03* Si_abun**3 + 4.38607436e+02 *Si_abun**2 -9.33878515e+00 *Si_abun +  7.90628624e+00
#    f2 = 0.0*Si_abun

### Biogenic silica parameterization for Fig. 6 in main text
#    f1 =  1028.61184527 *Si_abun**6 -2290.31344993 *Si_abun**5 +  1944.60677108*Si_abun**4  -755.71962734* Si_abun**3 +  134.09471848 *Si_abun**2   -2.8274206 *Si_abun +  7.5015002
#    f2 = 0.0*Si_abun

### Biogenic silica parameterization for Fig. 7 in main text
#    f1 =  205.87090281*Si_abun**3 -126.46613902*Si_abun**2  + 37.62121989* Si_abun  +  7.40879628 #3.6166122*Si_abun**3 + 0.30710427*Si_abun**2 + 4.83648973*Si_abun + 1.15605322
#    f2 = 0.0*Si_abun

    ## Gradual transition between Phanerozoic and Precambrian parameterizations:
    if t0<=0.541e9:
        flux = f1*1e12
    elif (t0>0.541e9)and(t0<1.041e9):
        flux = 1e12* (f1*(1.041e9-t0)/5e8  + f2*(t0-0.541e9)/5e8)
    else:
        flux = f2*1e12
        
    if Si_abun<0:
        flux =0.0
    return flux


## Abiotic silica flux parmaterization 
def Si_aBio_flux(Si_abun,t0,Temperature):
    
    Si_sat = 1000*10**(-8.476 - 485.24/Temperature - 2.268e-6*Temperature**2.0 + 3.068*np.log10(Temperature))
    Si_abun = Si_abun*1000.0 # convert mMol/kg from mol/kg
    
    ### Choose coefficient to reproduce [Si] record:
    Coefficient = 50e12
    Coefficient = 100e12
    #Coefficient = 70e12
    flux = Coefficient*(Si_abun/Si_sat)**4.4 
   
    return flux
