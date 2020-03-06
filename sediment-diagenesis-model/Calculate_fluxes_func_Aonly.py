import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt


def fluxes(DB0,DC0,phi0,omega0,KB0,cf0,k_RW0,c00,DA0,SA0,plotting0,FB0,Btau0,Aexp0):
    
    global DB,DC,phi,omega,KB,cf,k_RW,c0,DA,SA,plotting,FB,Btau,Aexp
    DB=DB0
    DC=DC0
    phi=phi0
    omega=omega0
    KB=KB0
    cf=cf0
    k_RW=k_RW0
    c0=c00
    DA=DA0
    SA=SA0
    plotting=plotting0
    FB=FB0
    Btau=Btau0
    Aexp = Aexp0
      
    x = np.linspace(0, 100, 1000) #default 1000

    y_a = 0.5*np.ones((4, x.size))
    y_a[2,:] = 500.0*y_a[2,:]
    
    y_b = np.zeros((6, x.size))

    res_a = solve_bvp(fun1, bc1, x, y_a,max_nodes=1000000,verbose=1) #100000 standard, 1000000 for when bumping SA to 0.9
    print (np.median(res_a.rms_residuals))
    
    x_plot = np.linspace(0, 100, 100000)
    y_plot_0 = res_a.sol(x_plot)[0]
    y_plot_1 = res_a.sol(x_plot)[1]
    y_plot_2 = res_a.sol(x_plot)[2]
    y_plot_3 = res_a.sol(x_plot)[3]

    if res_a.success == True:
        print (res_a.success)
    else:
        print( " ")
        print ("THIS ONE HAS A PROBLEM")
        print (" ")
    print ("Input parameters")
    print (DB0,DC0,phi0,omega0,KB0,cf0,k_RW0,c00,DA0,SA0,plotting0,FB0,Btau0,Aexp0)
    
    print ("Implied dissolved silica flux from/to ocean (+ source, - sink)")
    Dissolved_silica_back = - phi*y_plot_1[0]*DC
    print (Dissolved_silica_back)

    Advection_silica_sink = phi*omega*(y_plot_0[np.size(x_plot)-1]-y_plot_0[0]) ## loss at bottom - gain at top//
    #print(Advection_silica_sink)

    print ("A-silicate advection sink")
    Advection_A_sink = omega*(y_plot_2[np.size(x_plot)-1]-y_plot_2[0]) ## loss at bottom - gain at top//
    print(Advection_A_sink,omega*y_plot_2[np.size(x_plot)-1],omega*y_plot_2[0])

    print ("Implied A-silicate from/to ocean (+ source, - sink)") ## Problem because clays don't diffuse back into ocean! Or really diffuse at all!
    Dissolved_Asilica_back = - y_plot_3[0]*DA
    print (Dissolved_Asilica_back)
    print ("surface gains - advection losses",Dissolved_silica_back+Dissolved_Asilica_back - (Advection_A_sink+Advection_silica_sink))
    
    #### Sink partitioning 
    print ("Si loss via RW bottom of sediment column")
    OUT1=omega*(y_plot_2[np.size(x_plot)-1]) # cm/s * (umol/cm3 bulk) = umol/cm2/s
    print(OUT1)
    print("Si loss via dissolved Si bottom of sediment column")
    OUT2=phi*omega*(y_plot_0[np.size(x_plot)-1]) # cm3 liq/cm3 bulk * cm/s * (umol/cm3 liq) = umol/cm2/s
    print(OUT2)
    print ("Si loss via B Si bottom of sediment column")  # cm/s * (umol/cm3 bulk) = umol/cm2/s 
    OUT3=OUT1*0.0
    print(OUT3)
    OUT=OUT1+OUT2+OUT3
    print("RW sink,SiO2 sink, Bsink")
    print(OUT1/OUT,OUT2/OUT,OUT3/OUT)

    units = "a" # in umol/cm3
    #units = "b" # convert to weight % for plotting
    
    if plotting=="y":
        ### plotting

        if units =="b":
            y_plot_2 = 100*y_plot_2*1e-6*60/2.5
            y_plot_3 = 100*y_plot_3*1e-6*60/2.5

            
        plt.figure()
        Title = str(res_a.success)
        Title2 = res_a.message
        plt.title(Title+" "+Title2+" "+str(c0))

        plt.subplot(2,3,1)
        plt.title(Title+"                      ")
        plt.plot( y_plot_0,x_plot, label='c')
        plt.gca().invert_yaxis()
        #plt.legend()
        plt.ylabel("Depth (cm)")
        plt.xlabel('c (mM or umol/cm3 liq)')
            
        plt.subplot(2,3,2)
        plt.title(Title2)
        plt.plot( y_plot_1,x_plot, label='dc/dz')
        plt.gca().invert_yaxis()
        #plt.legend()
        plt.ylabel("Depth (cm)")
        plt.xlabel('dc/dz')

        plt.subplot(2,3,3)
        plt.title("                      "+str(c0))
        plt.plot( y_plot_2*0.0,x_plot, label='B')
        plt.gca().invert_yaxis()
        plt.ylabel("Depth (cm)")
        plt.xlabel('B (umol/cm3 sediment)')
        if units =="b":
            plt.xlabel('B (wt %)')
        
        plt.subplot(2,3,4)
        plt.plot( y_plot_3*0.0,x_plot, label='dB/dz')
        plt.gca().invert_yaxis()
        plt.ylabel("Depth (cm)")
        plt.xlabel('dB/dz')
        
        plt.subplot(2,3,5)
        plt.plot( y_plot_2,x_plot, label='A')
        plt.gca().invert_yaxis()
        plt.ylabel("Depth (cm)")
        plt.xlabel('A (umol/cm3 sediment)')
        if units =="b":
            plt.xlabel('A (wt %)')
            
        plt.subplot(2,3,6)
        plt.plot( y_plot_3,x_plot, label='dA/dz')
        plt.gca().invert_yaxis()
        plt.ylabel("Depth (cm)")
        plt.xlabel('dA/dz')

    return OUT1,OUT2,OUT3,OUT,y_plot_0,y_plot_2*0.0,y_plot_2
    
        
def fun1(x, y): 

    RW = 0.0*x
    RW =k_RW*( np.exp(-Aexp*x)+0.012/110.0) * ( y[0]/SA - 1) 
    RW[np.where(y[0]<SA)]=0.0 
    cofactor = np.exp( -Btau*x)  

    y[0][np.where(y[0]<0)] = 0
    y[2][np.where(y[2]<0)] = 0

    return np.vstack((y[1], RW/(DC*phi)+(omega/DC)*y[1],y[3],(omega/DA)*y[3]-RW/DA))

def bc1(ya, yb): 
    return np.array([ya[0]-c0, yb[1],ya[2],yb[3]]) 


            
