## This script loads the outputs from the sediment diagenesis models, and fits them 
## with simple functions (usually polynomials) for use in the coupled carbon-silica cycle model
## Specifically, these are used in the script Si_cycle_functions.py

import scipy.stats
import numpy as np
from Calculate_fluxes_func import fluxes
import matplotlib.pyplot as plt

### Choose whether fitting biogenic or non-biogenic case
[c0_range,RW_sink_1,RW_sink_2,RW_sink_3] = np.load("RW_array_no_Bflux.npy")
#[c0_range,RW_sink_1,RW_sink_2,RW_sink_3] = np.load("RW_array_Bflux.npy")
RW_total = RW_sink_1+RW_sink_2+RW_sink_3

[c0_rangeB,BSi_sink_1,BSi_sink_2,BSi_sink_3] = np.load("BSi_array_no_Bflux.npy")
#[c0_rangeB,BSi_sink_1,BSi_sink_2,BSi_sink_3] = np.load("BSi_array_Bflux.npy")
BSi_total = BSi_sink_1+BSi_sink_2+BSi_sink_3

### Choose appropriate polynomial fit (this needs to be done heuristically):
p = np.polyfit (c0_rangeB[0:15],BSi_total[0:15],3)
#p = np.polyfit (c0_rangeB[0:16],BSi_total[0:16],4)
#p = np.polyfit (c0_rangeB[0:18],BSi_total[0:18],6) #B
new_c0B = np.linspace(0.001,1.0,100)
#new_BSi = p[0]*new_c0B**3+p[1]*new_c0B**2+p[2]*new_c0B+p[3] #204.57442566 -140.36961165   38.55956195    4.47396617]
#new_BSi = p[0]*new_c0B**4+p[1]*new_c0B**3+p[2]*new_c0B**2+p[3]*new_c0B+p[4]
new_BSi = p[0]*new_c0B**3+p[1]*new_c0B**2+p[2]*new_c0B+p[3]
#new_BSi= p[0]*new_c0B**6+p[1]*new_c0B**5+p[2]*new_c0B**4+p[3]*new_c0B**3+p[4]*new_c0B**2+p[5]*new_c0B+p[6]

#q = np.polyfit (c0_range[0:26],RW_total[0:26],5) #B
#q = np.polyfit (c0_range[0:16],RW_total[0:16],6) #B
#q = np.polyfit (c0_range[0:16],RW_total[0:16],7) #B
#q = np.polyfit (c0_range[0:16],RW_total[0:16],4) #B
q = np.polyfit (c0_range[0:28],RW_total[0:28],3)
q = np.polyfit (c0_range[0:29],RW_total[0:29],3)
#q = np.polyfit (c0_range[0:15],RW_total[0:15],2) #B
new_c0 = np.linspace(0.001,2.7,100) #B
#new_c0 = np.linspace(0.01,1.0,50) #B
#q = np.polyfit (c0_range,RW_total,4) #noB
#new_c0 = np.linspace(0.02,1.5,50) #noB
#new_RW = q[0]*new_c0**4+q[1]*new_c0**3+q[2]*new_c0**2+q[3]*new_c0+q[4] #126.0539942  -211.47360998  103.27129996   -6.1339141     5.37447957
#new_RW = q[0]*new_c0**5+q[1]*new_c0**4+q[2]*new_c0**3+q[3]*new_c0**2+q[4]*new_c0+q[5]
#new_RW = q[0]*new_c0**6+q[1]*new_c0**5+q[2]*new_c0**4+q[3]*new_c0**3+q[4]*new_c0**2+q[5]*new_c0+q[6]
#new_RW = q[0]*new_c0**7 + q[1]*new_c0**6 + q[2]*new_c0**5 + q[3]*new_c0**4 + q[4]*new_c0**3 + q[5]*new_c0**2 + q[6]*new_c0 + q[7]
new_RW = q[0]*new_c0**3+q[1]*new_c0**2+q[2]*new_c0+ q[3]
#new_RW = q[0]*new_c0**2+q[1]*new_c0+q[2]

### alternative
from scipy.optimize import curve_fit

def funcfit(x,a,b,c,d,e):
    return a*np.tanh(x*b-c)+d*x+e

#popt,pcov = curve_fit(funcfit,c0_range[0:17],RW_total[0:17])

#c0_range=c0_range[0:18]
#RW_sink_1=RW_sink_1[0:18]
#RW_sink_2=RW_sink_2[0:18]
#RW_sink_3=RW_sink_3[0:18]
#RW_sink_4=RW_sink_4[0:18]
#RW_total=RW_total[0:18]

#BSi_sink_1=BSi_sink_1[0:18]
#BSi_sink_2=BSi_sink_2[0:18]
#BSi_sink_3=BSi_sink_3[0:18]
#BSi_sink_4=BSi_sink_4[0:18]
#BSi_total=BSi_total[0:18]

#print (popt)
plt.figure()
plt.plot(c0_range,RW_sink_1,'r',label=" proximal coast")
plt.plot(c0_range,RW_sink_2,'g',label="distal coast")
plt.plot(c0_range,RW_sink_3,'b',label="open ocean")
plt.plot(c0_range,RW_total,'k',label="total")
plt.plot(new_c0,new_RW,label="polynomial fit")
#plt.plot(new_c0,funcfit(new_c0,*popt),label="tanh fit")
plt.legend()
plt.ylabel('RW_sink')

plt.figure()
plt.plot(c0_rangeB,BSi_sink_1,'r',label="proximal coast")
plt.plot(c0_rangeB,BSi_sink_2,'g',label="distal coast")
plt.plot(c0_rangeB,BSi_sink_3,'b',label="open ocean")
plt.plot(c0_rangeB,BSi_total,'k',label="total")
plt.plot(new_c0B,new_BSi,label="polynomial fit")
plt.legend()
plt.ylabel('BSi_sink')
plt.show()

#plt.figure()
#plt.plot(c0_range,RW_sink_1,'r',label="tropical proximal coast")
#plt.plot(c0_range,RW_sink_2,'g',label="temperate proximal coast")
#plt.plot(c0_range,RW_sink_3,'b',label="distal coast")
#plt.plot(c0_range,RW_sink_4,'c',label="open ocean")
#plt.plot(c0_range,RW_total,'k',label="total")
#plt.plot(new_c0,new_RW,label="polynomial fit")
#plt.legend()
#plt.ylabel('RW_sink')

#plt.figure()
#plt.plot(c0_rangeB,BSi_sink_1,'r',label="tropical proximal coast")
#plt.plot(c0_rangeB,BSi_sink_2,'g',label="temperate proximal coast")
#plt.plot(c0_rangeB,BSi_sink_3,'b',label="distal coast")
#plt.plot(c0_rangeB,BSi_sink_4,'c',label="open ocean")
#plt.plot(c0_rangeB,BSi_total,'k',label="total")
#plt.plot(new_c0,new_BSi,label="polynomial fit")
#plt.legend()
#plt.ylabel('BSi_sink')
#plt.show()
