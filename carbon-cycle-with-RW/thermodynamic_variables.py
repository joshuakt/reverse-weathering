import pylab

# The function Sol_prod calculates carbonate solubility product as 
# a function of temperature. See Appendix A in Krissansen-Totton et
# al. (2018; PNAS) for further information. Expression originally 
# described in appendix G in Pilson 1998 # (reproduces table G.1 
# for salinity=35).

def Sol_prod(T): 
    bo=-.77712
    b1=0.0028426
    b2=178.34
    co=-.07711
    do=.0041249
    S=35
    (bo+b1*T+b2/T)*S**0.5+co*S+do*S**1.5
    logK0=-171.9065-0.077993*T+2839.319/T+71.595*pylab.log10(T) 
    logK=logK0+(bo+b1*T+b2/T)*S**0.5+co*S+do*S**1.5
    return 10**logK

# The function equil_cont calculates carbon chemistry equilibrium constants 
# as a function of temperature:
def equil_cont(T):
    pK1=17.788 - .073104 *T - .0051087*35 + 1.1463*10**-4*T**2
    pK2=20.919 - .064209 *T - .011887*35 + 8.7313*10**-5*T**2
    H_CO2=pylab.exp(9345.17/T - 167.8108 + 23.3585 * pylab.log(T)+(.023517 - 2.3656*10**-4*T+4.7036*10**-7*T**2)*35)
    return [pK1,pK2,H_CO2]