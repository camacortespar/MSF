
import numpy as np
import spiceypy as spice

import planetary_data as pd
import tools as _t
import spice_tools as st
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts


#central body
cb=pd.earth

#total time
tspan=3600*24*365.0

dt=1000

#central body
cb=pd.earth

date0='2020-02-23'

h=30.0e-3
w=35.0e-3
A=h*w

if __name__ == '__main__':

        #initial conditions
        state0=_t.tle2coes('HJ-2A.txt')

        state1=[42095.0,0.81818,28.5,180.0,298.2253,357.857]

        #perts dictionary
        perts=null_perts()

        #add solar radiation preasure
        perts['srp']=True
        perts['a_srp']=A
        perts['CR']=1.0
        mass0=7000.0


        #create orbit propagator instance geostacionary satelite and ISS
        op0=OP(state0,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,mass0=mass0,propagator='dopri5')
        op1=OP(state1,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,mass0=mass0,propagator='dopri5')

        #op0.calculate_coes(parallel=False)
        #op0.plot_coes(days=True,show_plot=True,rel=True,title='Orbita HJ_2A con SRP')

        op1.calculate_coes(parallel=False)
        op1.plot_coes(days=True,show_plot=True,rel=True,title='Orbita satélite Molniya con SRP')

        _t.plot_n_orbits([op0.rs,op1.rs],labels=['HJ_2A','Molniya Orbit'],show_plot=True,title='Orbita satélite Molniya y HJ-2A con SRP')