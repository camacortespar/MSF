
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
tspan=60623.4

#time step
dt=1

date0='5/12/2020-10:19:19'

mass0=1457

if __name__ == '__main__':

        #initial conditions
        iss_coes=_t.tle2coes('NOAA15.txt',degres=True)

        #state0=[42164.0,0.001,0.0,0.0,0.0,0.0]

        #perts dictionary
        perts=null_perts()
        perts['J2']=True

        perts['aero']=True
        perts['Cd']=2.0
        perts['A']=0.00031 #km²

        perts['srp']=True
        perts['A_srp']=0.00031
        perts['CR']=1

        #perts['BSTAR']=True
        #perts['B*']=0.000020351

        #add n_body gravity perturbation
        perts['n_bodies']=[pd.moon]+[pd.sun]+[pd.jupiter]+[pd.venus]+[pd.mars]+[pd.mercury]+[pd.saturn]+[pd.neptune]+[pd.uranus]

       
        op_ISS=OP(iss_coes,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,mass0=mass0,propagator='lsoda')

        op_ISS.calculate_coes(parallel=False)

        op_ISS.printcoes()
        #op_ISS.plot_coes(days=True,show_plot=False,title='Coes ISS perturbación de la luna')

        #rs_m=np.array(op_ISS.perts['n_bodies'][0]['states'][:,:3])
        #rs_m2=np.array(op_ISS.perts['n_bodies'][1]['states'][:,:3])

        #op_ISS.plot_3d(show_plot=True)

       # _t.plot_n_orbits([op_ISS.rs],labels=['ISS'],show_plot=True)
