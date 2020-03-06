## This script uses a MCMC approach to find sediment diagenesis parameters to fit modern conditions.
## Code needs to be run separately for proximal, distal, and pelagic conditions
## Also, two cases for T&D modern RW flux and Rahman et al. modern RW flux.

import emcee
import numpy
import matplotlib.pyplot as pl
from corner import corner
import corner
import numpy as np
from Calculate_fluxes_func import fluxes 

##########################################################
########## Case 1: T&D Initial Conditions ################
##########################################################

## Fit for total proximal coastal
obs_RW_flux = 1.0 #Tmol/yr
er_RW_flux = 0.25 # Tmol/yr
obs_B_flux = 3.0 # Tmol/yr
er_B_flux = 0.25 # Tmol/yr
obs_cDEEP = 0.5 # mM
er_cDEEP = 0.4 # mM
obs_BDEEP = 500 # um/bulk cm3
er_BDEEP = 1000 # um/bulk cm3

DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
phi=0.90 ## porosity constant
omega_1 =6e-9 #weighted averae of tropical and temparate
cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
c0=0.050 # dissolved silica in ocean
DA = 1e-15#5e-8
area_fraction = 0.0052
###################################

## Fit for distal coastal box
###########################
#obs_RW_flux = 0.5 #Tmol/yr
#er_RW_flux = 0.25 # Tmol/yr
#obs_B_flux = 3.0 # Tmol/yr
#er_B_flux = 0.25 # Tmol/yr
#obs_cDEEP = 0.5 # mM
#er_cDEEP = 1.0 # mM
#obs_BDEEP = 500 # um/bulk cm3
#er_BDEEP = 1000 # um/bulk cm3

#DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
#DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
#phi=0.90 ## porosity constant
#omega_1 =9.8e-10 #midway 3.1e-9, 3.1e-10 in log space
#cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
#c0=0.050 # dissolved silica in ocean
#DA = 1e-15 #DB
#area_fraction = 0.076
###############################

## Fit for pelagic (deep ocean)
###########################
obs_RW_flux = 0.0 #Tmol/yr
er_RW_flux = 0.25 # Tmol/yr
obs_B_flux = 3.0 # Tmol/yr
er_B_flux = 0.05#0.25 # Tmol/yr
obs_cDEEP = 0.5 # mM
er_cDEEP = 2.0 # mM
obs_BDEEP = 500 # um/bulk cm3
er_BDEEP = 2000 # um/bulk cm3

DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
phi=0.90 ## porosity constant
omega_1 =1e-12  #global median good enough?
cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
c0=0.050 # dissolved silica in ocean
DA = 1e-15#DB
area_fraction = 0.919
###############################

#############################################################
########## Case 2: Rahman Initial Conditions ################
#############################################################

## Fit for total proximal coastal
#obs_RW_flux = 4.5 #Tmol/yr
#er_RW_flux = 0.25 # Tmol/yr
#obs_B_flux = 3.7 # Tmol/yr
#er_B_flux = 0.25 # Tmol/yr
#obs_cDEEP = 0.5 # mM
#er_cDEEP = 0.4 # mM
#obs_BDEEP = 500 # um/bulk cm3
#er_BDEEP = 1000 # um/bulk cm3

#DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
#DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
#phi=0.90 ## porosity constant
#omega_1 =6e-9 #weighted averae of tropical and temparate
#cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
#c0=0.050 # dissolved silica in ocean
#DA = 1e-15#5e-8
#area_fraction = 0.0052
###################################

#### fit for distal coastal
#obs_RW_flux = 0.01 #Tmol/yr
#er_RW_flux = 0.5 # Tmol/yr
#obs_B_flux = 3.0 # Tmol/yr
#er_B_flux = 0.5 # Tmol/yr
#obs_cDEEP = 0.5 # mM
#er_cDEEP = 1.0 # mM
#obs_BDEEP = 500 # um/bulk cm3
#er_BDEEP = 1000 # um/bulk cm3

#DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
#DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
#phi=0.90 ## porosity constant
#omega_1 =9.8e-10 #midway 3.1e-9, 3.1e-10 in log space
#cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
#c0=0.050 # dissolved silica in ocean
#DA = 1e-15 #DB
#area_fraction = 0.076
########################################

# deep pelagic ocean deposition
#obs_RW_flux = 0.0 #Tmol/yr
#er_RW_flux = 0.25 # Tmol/yr
#obs_B_flux = 1.0 # Tmol/yr
#er_B_flux = 0.2 # Tmol/yr
#obs_cDEEP = 0.5 # mM
#er_cDEEP = 2.0 # mM
#obs_BDEEP = 500 # um/bulk cm3
#er_BDEEP = 2000 # um/bulk cm3

#DB=5e-8 ## effective diffusion coefficient for biogenic silica, non-zero because bioturbation (Schink numbers)
#DC=4e-6 ## diffusion coefficient accounting for porosity and tortuoisty for dissolved silica (Schink numbers)
#phi=0.90 ## porosity constant
#omega_1 =1e-12  #global median good enough?
#cf=1.0 ## saturation state for biogenic silica (from literature e.g. Rickert 2002)
#c0=0.050 # dissolved silica in ocean
#DA = 1e-15#DB
#area_fraction = 0.919
###############################

scale = (365*24*60*60) *(0.71*5.1e18)/1e6# uMol/cm2/s

def LnLike(x):
  
  if (x[0] > -8.0 ) or (x[0]< -12.0): ### logo10( biogenic dissolution rate (Schink))
    return -numpy.inf    
  if (x[1] > -5.4 ) or (x[1]< -10.0): ##log10(k_RW, A reaction rate)
    return -numpy.inf    
  if (x[2] > 1.0 ) or (x[2] < -3.0): ## log B_Tau
    return -numpy.inf #added
  if (x[3] > -5.0) or (x[3] < -10.0): #log10( FB_change)
    return -numpy.inf #added
  if (x[4] > 1.0) or (x[4] < 0.2): #SA
    return -numpy.inf #added
  if (x[5] > 0) or (x[5] < -3): #logo10(ATau)
    return -numpy.inf #added

  try:
      [RW_sink,dSI_sink,BSi_sink,Total_sink,c_profile,B_profile,A_profile] = fluxes(DB,DC,phi,omega_1,10**(x[0]),cf,10**(x[1]),c0,DA,x[4],"n",10**(x[3]),10**(x[2]),10**(x[5])) 
  except:
      return -numpy.inf

  if np.any(c_profile)>1.1:
      return -numpy.inf

  model_RW_flux =area_fraction*np.array(RW_sink)*scale/1e12 ## Tropical box
  model_B_flux =numpy.max([area_fraction*np.array(BSi_sink)*scale/1e12,1e-20]) ## max needed to ensure -1e-20 values don'y mess up log
  model_cDEEP = c_profile[90000]
  model_BDEEP = B_profile[90000]

  ## Additional ratio constraints
  Bflux_model = area_fraction*(10**(x[3]))*scale/1e12
  B_advection_source = area_fraction*B_profile[0]*omega_1*scale/1e12

  try:
      model_Bburial_to_Bflux = numpy.log10(model_B_flux / (Bflux_model+B_advection_source))
  except:
      return -numpy.inf

  sediment_flux=omega_1*2.5 # cm/s * g/cm3 = g sed/cm2/s
  SiO2_flux = (10**(x[3]) + B_profile[0]*omega_1)*1e-6*60 # umol/cm2/s * 10-6 mol/umol * 60 g/mol =g SiO/cm2/s, um/cm3 * cm/s * mol/umol * g/mol = g SiO2/cm2/2
  if sediment_flux <SiO2_flux:
      return -numpy.inf

  obs_ratio = -1.0 # ensures about 10% of gross flux is permanently buried
  er_obs_ratio = 0.3

  log_terms=numpy.sum(numpy.log(2*numpy.pi*er_RW_flux**2))   + numpy.sum(numpy.log(2*numpy.pi*er_B_flux**2)) +numpy.sum(numpy.log(2*numpy.pi*er_cDEEP**2)) + numpy.sum(numpy.log(2*numpy.pi*er_BDEEP**2)) + numpy.sum(numpy.log(2*numpy.pi*er_obs_ratio**2))
  log_like =  -0.5 * (  numpy.sum((obs_RW_flux-model_RW_flux)**2/er_RW_flux**2) + numpy.sum((obs_B_flux-model_B_flux)**2/er_B_flux**2) + numpy.sum((obs_cDEEP-model_cDEEP)**2/er_cDEEP**2) + numpy.sum((obs_BDEEP-model_BDEEP)**2/er_BDEEP**2)+ numpy.sum((model_Bburial_to_Bflux-obs_ratio)**2/er_obs_ratio**2) )

  return log_like

# Set up the sampler
nwalk = 50
ndim = 6 
nsteps = 500

p0 = numpy.vstack([[-12+4.0*numpy.random.random() for i in range(nwalk)],
                [-8.0+6.0*numpy.random.random() for i in range(nwalk)],
                [-3.0+4.0*numpy.random.random() for i in range(nwalk)],
                [-10+5.0*numpy.random.random() for i in range(nwalk)],
                [0.2+0.8*numpy.random.random() for i in range(nwalk)],
                [-3.0+3.0*numpy.random.random() for i in range(nwalk)]]).T # added
sampler = emcee.EnsembleSampler(nwalk, ndim, LnLike,threads=30) #8 safer maybe
                                #args = (t, y, error), 
                                #pool = emcee.interruptible_pool.InterruptiblePool()) #for parallel
pos, lnprob, rstate=sampler.run_mcmc(p0, nsteps)

numpy.save('chain_test',sampler.chain[:,:,:]) #originally chain
numpy.save('lnprob_test',sampler.lnprobability) #originally lnprob

chain=numpy.load('chain_test.npy') #originally chain
lnprob=numpy.load('lnprob_test.npy') #originally lnprob

# Plot the chains
fig, ax = pl.subplots(6)
for n in range(nwalk):
  ax[0].plot(chain[n,:,0])
  ax[1].plot(chain[n,:,1])
  ax[2].plot(chain[n,:,2]) 
  ax[3].plot(chain[n,:,3])
  ax[4].plot(chain[n,:,3])
  ax[5].plot(chain[n,:,3])

#find highest likelihood run
logprob=numpy.array(lnprob)
values=chain
ii,jj = numpy.unravel_index(logprob.argmax(), logprob.shape)
print ("indeces for best",ii,jj)
print ("loglikelihood and values",logprob[ii,jj],values[ii,jj,:])
print ("Check",LnLike(values[ii,jj,:]))

# Plot the corner plot, discarding the first 100 steps as burn-in
production = chain[:,0:,:]
s = production.shape
flatchain = production.reshape(s[0] * s[1], s[2])
#flatchain=sampler.flatchain

# more limited corner plots
corner.corner(flatchain[:,[0,1,2,3,4,5]], quantiles=[0.16, 0.5, 0.84],labels=["KB", "K_RW", "Btau","FB","SA","ATau"])

flatchain2=flatchain
from matplotlib import rc

ab, bc, cd,de,ef,fg= map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*numpy.percentile(flatchain, [16, 50, 84],axis=0)))
print ("median values with errors", numpy.array([ab, bc, cd,de,ef,fg]))
#print("Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction)))
print ("confidence intervals")
map(lambda v: (v[0], v[1], v[2]),zip(*numpy.percentile(flatchain, [5, 50, 95],axis=0)))

[RW_sink_best,dSI_sink_best,BSi_sink_best,Total_sink_best,c_profile_best,B_profile_best,A_profile_best] = fluxes(DB,DC,phi,omega_1,10**(values[ii,jj,0]),cf,10**(values[ii,jj,1]),c0,DA,values[ii,jj,4],"y",10**(values[ii,jj,3]),10**(values[ii,jj,2]),10**(values[ii,jj,5]))
print ("KB",10**(values[ii,jj,0]),"k_RW",10**(values[ii,jj,1]),"BTau",10**(values[ii,jj,2]),"FB_1",10**(values[ii,jj,3]),"SA",(values[ii,jj,4]),"ATau",10**(values[ii,jj,5]))

RW_sink_best = RW_sink_best*area_fraction*scale/1e12
dSI_sink_best = dSI_sink_best*area_fraction*scale/1e12
BSi_sink_best = BSi_sink_best*area_fraction*scale/1e12
Total_sink_best = Total_sink_best*area_fraction*scale/1e12
print (RW_sink_best,dSI_sink_best,BSi_sink_best,Total_sink_best)

import pylab
legend_counter=0
for x_ex in flatchain[numpy.random.randint(len(flatchain), size=10)]:
    [a,b,c,d,e,f,g]=fluxes(DB,DC,phi,omega_1,10**(x_ex[0]),cf,10**(x_ex[1]),c0,DA,x_ex[4],"y",10**(x_ex[3]),10**(x_ex[2]),10**(x_ex[5]))
    legend_counter=legend_counter+1

mega_output=[]
for x_ex in flatchain[numpy.random.randint(len(flatchain), size=40)]:
    [aa0,bb0,cc0,dd0,ee0,ff0,gg0]=fluxes(DB,DC,phi,omega_1,10**(x_ex[0]),cf,10**(x_ex[1]),c0,DA,x_ex[4],"n",10**(x_ex[3]),10**(x_ex[2]),10**(x_ex[5]))
    outputs = [ee0,ff0,gg0]
    mega_output.append(outputs)

mega_output=numpy.array(mega_output)

import scipy.stats
from scipy.interpolate import interp1d
confidence_c=scipy.stats.scoreatpercentile(mega_output[:,0,:],[5.0,50,95], interpolation_method='fraction',axis=0)
confidence_B=scipy.stats.scoreatpercentile(mega_output[:,1,:],[5.0,50,95], interpolation_method='fraction',axis=0)
confidence_A=scipy.stats.scoreatpercentile(mega_output[:,2,:],[5.0,50,95], interpolation_method='fraction',axis=0)

x_plot = np.linspace(0, 100, 100000)
import pylab as plt

plt.figure()
plt.subplot(1,3,1)
plt.plot(c_profile_best, x_plot, label='c')
pylab.fill_betweenx(x_plot, confidence_c[0], confidence_c[2], where=confidence_c[2]>= confidence_c[0], facecolor='red', alpha='0.4')
plt.gca().invert_yaxis()
plt.ylabel("Depth (cm)")
plt.xlabel('c (mM or umol/cm3 liq)')

plt.subplot(1,3,2)
plt.plot(B_profile_best, x_plot, label='B')
pylab.fill_betweenx(x_plot, confidence_B[0], confidence_B[2], where=confidence_B[2]>= confidence_B[0],facecolor='red', alpha='0.4')
plt.gca().invert_yaxis()
plt.ylabel("Depth (cm)")
plt.xlabel('B (umol/cm3 sediment)')

plt.subplot(1,3,3)
plt.plot(A_profile_best, x_plot, label='A')
pylab.fill_betweenx(x_plot, confidence_A[0], confidence_A[2], where=confidence_A[2]>= confidence_A[0], facecolor='red', alpha='0.4')
plt.gca().invert_yaxis()
plt.ylabel("Depth (cm)")
plt.xlabel('A (umol/cm3 sediment)')

[RW_sink_best,dSI_sink_best,BSi_sink_best,Total_sink_best,c_profile_best,B_profile_best,A_profile_best] = fluxes(DB,DC,phi,omega_1,10**(values[ii,jj,0]),cf,10**(values[ii,jj,1]),c0,DA,values[ii,jj,4],"y",10**(values[ii,jj,3]),10**(values[ii,jj,2]),10**(values[ii,jj,5]))
print ("KB",10**(values[ii,jj,0]),"k_RW",10**(values[ii,jj,1]),"BTau",10**(values[ii,jj,2]),"FB_1",10**(values[ii,jj,3]),"SA",(values[ii,jj,4]),"ATau",10**(values[ii,jj,5]))
RW_sink_best = RW_sink_best*area_fraction*scale/1e12
dSI_sink_best = dSI_sink_best*area_fraction*scale/1e12
BSi_sink_best = BSi_sink_best*area_fraction*scale/1e12
Total_sink_best = Total_sink_best*area_fraction*scale/1e12
print (RW_sink_best,dSI_sink_best,BSi_sink_best,Total_sink_best)

pl.show()