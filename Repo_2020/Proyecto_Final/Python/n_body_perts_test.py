
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
tspan=3600*24*100.0

#time step
dt=1000

#central body
cb=pd.earth

date0='2020-02-23'

if __name__ == '__main__':

        #initial conditions
        iss_coes=_t.tle2coes('iss.txt')

        state0=[42164.0,0.001,0.0,0.0,0.0,0.0]

        #perts dictionary
        perts=null_perts()

        #add lunar gravity perturbation
        perts['n_bodies']=[pd.moon]

        #create orbit propagator instance geostacionary satelite and ISS
        op0=OP(state0,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,propagator='dopri5')
        op_ISS=OP(iss_coes,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,propagator='dopri5')

        op_ISS.calculate_coes(parallel=False)
        op_ISS.plot_coes(days=True,show_plot=False,title='Coes ISS perturbación de la luna')

        op0.calculate_coes(parallel=False)
        op0.plot_coes(days=True,show_plot=False,title='Coes geostacionario perturbación de la luna')

        rs_m=np.array(op0.perts['n_bodies'][0]['states'][:,:3],dtype=object)

        _t.plot_n_orbits([op_ISS.rs,op0.rs,rs_m],labels=['ISS','GEO','Moon'],show_plot=True)