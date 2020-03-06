# This function runs the sediment diagenesis model for a range of ocean boundary conditions
# i.e. different dissolved silica concentrations. This can be done for all the scenarios
# in the main text (Fig. 3 - 7). Additionally, the sediment diagenesis model can be run
# either with or without biogenic silica precipitation, representing cases in the Phanerozoic
# and Precambrian, respectively.

import scipy.stats
import numpy as np
import matplotlib.pyplot as plt

############################################################################################
####### OPTIONS ############################################################################
which_case = "y" # "n" for Precambrian, no biogenic silica precipitation, "y" Phanerozoic cases with biogenic silica precipitation
which_fig = 3 # options 3,4,5,6,7 for reproducing silica functions for Fig. 3-7 in main text
############################################################################################

if which_case == "n":
    from Calculate_fluxes_func_Aonly import fluxes
else:
    from Calculate_fluxes_func import fluxes

######
#distal coastal parameters
RW_sink_1=[]
dSI_sink_1=[]
BSi_sink_1=[]
Total_sink_1=[]
omega_1 =6e-09 #sed rate
area_fraction1 = 0.0052

##########
#proximal coastal parameters
RW_sink_2=[]
dSI_sink_2=[]
BSi_sink_2=[]
Total_sink_2=[]
omega_2 = 9.8e-10  #sed rate
area_fraction2 = 0.076

###########
#distal coastal parameters
RW_sink_3=[]
dSI_sink_3=[]
BSi_sink_3=[]
Total_sink_3=[]
omega_3 = 1e-12 #sed rate
area_fraction3 = 0.919
###########

## biogenic flux boundary condition
if which_fig == 3:
    FB_1 = 2.63485056405e-06 #Fig. 3 parameterization
    FB_2= 1.24041815806e-06 ##Fig. 3 parameterization
    FB_3 =  2.89332859848e-08 #Fig. 3 parameterization
elif which_fig == 4:
    FB_1 = 9.95588576597e-06 # Fig. 4 parameterization
    FB_2 =2.95681747708e-06 ## Fig. 4 parameterization
    FB_3 = 3.59682053804e-08 # Fig. 4 parameterization
elif (which_fig == 5)or(which_fig ==6):
    FB_1 = 1.10653408373e-06 # Fig. 5, 6 parameterization
    FB_2 =2.30004791083e-06 ## Fig. 5, 6 parameterization
    FB_3 = 1.44600620805e-08 # Fig. 5, 6 parameterization
else: 
    FB_1 = 5.2e-6 # Fig. 7 parameterization
    FB_2 = 2.8e-7 # Fig. 7 parameterization
    FB_3 = 1.4673416466353787e-08 # Fig. 7 parameterization


if which_case == "n":
    c0_range=np.linspace(0.001,2.5,20) # Precambrian, no biogenic silica precipitation case
else:     
    c0_range=np.linspace(0.05,0.98,20) #for Phanerozoic, biogenic silica precipitation case

scale = (365*24*60*60) *(0.71*5.1e18)/1e6# uMol/cm2/s to mol/yr conversion

c_er_array = []
for i in range(0,np.size(c0_range)):
    print(i)
    c0=c0_range[i]
    
    ############
    if which_case == "n":
        FB_1= 1e-30
        FB_2 = 1e-30
        FB_3 = 1e-30
        FB_4 = 1e-30

    DA=1e-15  
    if which_fig ==3:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_1,2.08741263774e-09,1.0,1e-6,c0,DA,0.25,"y",FB_1,0.352416172635,0.7)  # Fig. 3 paramaterization
    elif which_fig ==4:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_1,7.64442360835e-09,1.0,3.86326148881e-06,c0,DA,0.25,"y",FB_1,0.108806868435,0.7) # Fig. 4 parameterization
    elif which_fig ==5:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes( 5e-08,4e-06,0.9,omega_1,4.50926951726e-10,1.0,9.81578769698e-09,c0,DA,0.229256976887,"y",FB_1,0.0316702594575,0.00497280027478) ## Fig. 5 parameterization
    elif which_fig ==6:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes( 5e-08,4e-06,0.9,omega_1,4.50926951726e-10,1.0,4e-6,c0,DA,0.3,"y",FB_1,0.316702594575,0.00497280027478) # Fig. 6 parameterization
    else:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_1, 8.130256979767599e-09,1.0,1.7852225077073948e-06,c0,DA,0.6226001281387643,"y",FB_1,0.468277365622558,0.12045681021203046) # Fig. 7 parameterization
    
    if (A_profile[99999]>22000)or(RW_sink>1e5):
        RW_sink = omega_1*21000+0*RW_sink
        Total_sink = RW_sink + dSI_sink + BSi_sink
    if (np.any(B_profile) > 50000)or(BSi_sink>1e5):
        BSi_sink = omega_1 * 50000 + BSi_sink *0.0
        Total_sink = RW_sink + dSI_sink + BSi_sink

    B_in_1 = FB_1 + omega_1*B_profile[0]
    print ("Bout/Bin",BSi_sink/B_in_1)
    RW_sink_1.append(RW_sink)
    dSI_sink_1.append(dSI_sink)
    BSi_sink_1.append(BSi_sink)
    Total_sink_1.append(Total_sink)
    
    if which_fig ==3:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_2,4.42057892626e-09,1.0,1e-7,c0 ,DA, 0.25, "y", FB_2, 1.07490672549,0.7) # Fig. 3 paramaterization
    elif which_fig ==4:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_2,7.10418443447e-09,1.0,3.00320675323e-09,c0,DA,0.25,"y", FB_2, 0.0149086868855,0.7) # Fig. 4 parameterization
    elif which_fig ==5:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_2,4.58460800854e-09,1.0,7.13366535715e-09,c0,DA,0.403137922379,"y",FB_2,0.00161597465734,0.162195505982) # Fig. 5 parameterization
    elif which_fig ==6:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_2,4.58460800854e-09,1.0,2e-6,c0,DA,0.3,"y",FB_2,0.00161597465734,0.162195505982) #Fig. 6 parameterization
    else:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_2,7.819474009452296e-09,1.0,5.404224878197497e-07,c0,DA,0.8928384136524511,"y",FB_2,0.37169331041072834,0.09609291652633287) # Fig. 7 parameterization

    if (A_profile[99999]>22000)or(RW_sink>1e5):
        RW_sink = omega_2*21000+0*RW_sink
        Total_sink = RW_sink + dSI_sink + BSi_sink
    if (np.any(B_profile) > 50000)or(BSi_sink>1e5):
        BSi_sink = omega_2 * 50000 + BSi_sink *0.0
        Total_sink = RW_sink + dSI_sink + BSi_sink

    B_in_2 = FB_2 + omega_2*B_profile[0]
    print ("Bout/Bin",BSi_sink/B_in_2)
    RW_sink_2.append(RW_sink)
    dSI_sink_2.append(dSI_sink)
    BSi_sink_2.append(BSi_sink)
    Total_sink_2.append(Total_sink)

    if which_fig ==3:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_3,3.86675081447e-12,1.0,1e-8,c0, DA,0.25,"y",FB_3,1.52066018874,0.7)  # Fig. 3 paramaterization
    elif which_fig ==4:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_3,3.39637963621e-11,1.0,2.99599076241e-07,c0,DA,0.25,"y",FB_3,4.59903882772,0.7) # Fig. 4 paramaterization
    elif which_fig ==5:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_3,1.32765334506e-12,1.0,1.81433970977e-07,c0,DA,0.628246936439,"y" ,FB_3,2.35741275621,0.00353631842707) # Fig. 5 parameterization
    elif which_fig ==6:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_3,1.32765334506e-12,1.0,1e-6,c0,DA,0.8,"y" ,FB_3,2.35741275621,0.00353631842707) # Fig. 6 parameterization
    else:
        RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile = fluxes(5e-08,4e-06,0.9,omega_3,2.0762179135497078e-12,1.0,5.7453143319964955e-09,c0,DA,0.820371634862408,"y",FB_3,3.6234270796809547,0.002380470071604743) # Fig. 7 parameterization
        
    if (A_profile[99999]>22000)or(RW_sink>1e5):
        print ("ERROR HERE WITH A_profile")
        c_er_array.append(c0)
        RW_sink = omega_3*21000+0*RW_sink
        Total_sink = RW_sink + dSI_sink + BSi_sink
    if (np.any(B_profile) > 50000)or(BSi_sink>1e5):
        BSi_sink = omega_3 * 50000 + BSi_sink *0.0
        Total_sink = RW_sink + dSI_sink + BSi_sink
    
    B_in_3 = FB_3 + omega_3*B_profile[0]
    print ("Bout/Bin",BSi_sink/B_in_3)
    RW_sink_3.append(RW_sink)
    dSI_sink_3.append(dSI_sink)
    BSi_sink_3.append(BSi_sink)
    Total_sink_3.append(Total_sink)


    Total_B_in = (area_fraction1*B_in_1 + area_fraction2*B_in_2 + area_fraction3*B_in_3 )*scale/1e12
    Total_B_out = (area_fraction1*BSi_sink_1[i] + area_fraction2*BSi_sink_2[i] + area_fraction3*BSi_sink_3[i] )*scale/1e12
    print ("Total_B_in",Total_B_in,"Total_B_out",Total_B_out,"Ratio",Total_B_out/Total_B_in)

print (c_er_array)

RW_sink_1=area_fraction1*np.array(RW_sink_1)*scale/1e12 # Tmol/yr
RW_sink_2=area_fraction2*np.array(RW_sink_2)*scale/1e12 # Tmol/yr   
RW_sink_3=area_fraction3*np.array(RW_sink_3)*scale/1e12 # Tmol/yr

BSi_sink_1=area_fraction1*np.array(BSi_sink_1)*scale/1e12 # Tmol/yr 
BSi_sink_2=area_fraction2*np.array(BSi_sink_2)*scale/1e12 # Tmol/yr 
BSi_sink_3=area_fraction3*np.array(BSi_sink_3)*scale/1e12 # Tmol/yr 

## Save outputs for making polynomial fits
RW_array=np.array([c0_range,RW_sink_1,RW_sink_2,RW_sink_3])
if which_case == "n":
    np.save("RW_array_no_Bflux",RW_array)
else:
    np.save("RW_array_Bflux.npy",RW_array)
    
BSi_array=np.array([c0_range,BSi_sink_1,BSi_sink_2,BSi_sink_3])
if which_case == "n":
    np.save("BSi_array_no_Bflux",BSi_array)
else:
    np.save("BSi_array_Bflux",BSi_array)
    
plt.figure()
plt.plot(c0_range,RW_sink_1,'r') 
plt.plot(c0_range,RW_sink_2,'g')     
plt.plot(c0_range,RW_sink_3,'b')
plt.plot(c0_range,RW_sink_1+RW_sink_2+RW_sink_3,'k')
plt.ylabel('RW_sink')

plt.figure()
plt.plot(c0_range,BSi_sink_1,'r') 
plt.plot(c0_range,BSi_sink_2,'g')     
plt.plot(c0_range,BSi_sink_3,'b')
plt.plot(c0_range,BSi_sink_1+BSi_sink_2+BSi_sink_3,'k')
plt.ylabel('BSi_sink')
plt.show()
