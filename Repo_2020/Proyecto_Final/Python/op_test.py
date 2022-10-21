import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

#for sys import path
#path.append(~/Escritorio/Simulacion/SPICE)
#from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd

cb=pd.sun


if __name__== '__main__':

    #initial conditions
    r_mag=cb['radius']+cb['radius']*4 #km
    v_mag=np.sqrt(cb['mu']/r_mag) #km/s
    
    #initial conditions and velocity vectors
    r0=[r_mag,r_mag*0.01,r_mag*-0.5]
    v0=[0,v_mag,v_mag*0.8]
    
    #time span
    tspan=6*3600*24.0
    
    #time step
    dt=100.0

    op=OP(r0,v0,tspan,dt,cb=cb)
    op.propagate_orbit()
    op.plot_3d(show_plot=True)
    
