import numpy as np
import datetime
import spiceypy as spice
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
import tools as _t
import spice_tools as st
import lambert_tools as lt
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

dt=1000.0

#central body
cb=pd.sun

#initial conditions 
date0='2016-04-25'
datef='2024-08-15'

#reference frame
FRAME='ECLIPJ2000'

#center of reference frame
OBSEVER='SUN'

#header of csv files 
HEADER='t,rx,ry,rz'

if __name__ == '__main__':
	
	#load metakernel for solar system ephemerides
	spice.furnsh('spice_data/solar_system_kernel.mk')

	#convert start and end dates to seconds past J2000
	et0=spice.utc2et(date0)
	etf=spice.utc2et(datef)

	#calculate transfer time 
	transfer_time=etf-et0

	#time array for earth and venus
	tiem_arr=np.linspace(et0,etf,10000)

	#calculate earth and venus state vector at initial time
	states_earth=st.get_ephemerides_data('JUPITER BARYCENTER',tiem_arr,FRAME,OBSEVER)
	states_venus=st.get_ephemerides_data('SATURN BARYCENTER',tiem_arr,FRAME,OBSEVER)

	#space craft initial position vector
	r0=states_earth[0,:3]

	#calculate venus position vector at final time
	rf=states_venus[-1,:3]

	#calculate spacecraft velocity vectors via lamberrts solution
	v0,vf=lt.lamberts_problem_universal_variables(r0,rf,transfer_time,mu=cb['mu'])

	#initial state vector for spacecraft
	state0_sc=r0.tolist()+v0.tolist()

	#propagate spacecraft orbit
	op_sc=OP(state0_sc,transfer_time,dt,cb=cb)

	#
	#
	#
	lt.print_delta_v(v0,vf,states_earth[1,3:],states_venus[-1,3:])

	_t.plot_n_orbits([states_earth[:,:3],states_venus[:,:3],op_sc.rs]
		,labels=['Earth','Mars','Spacecraft']
		,cb=cb,show_plot=True,title='Jupiter to Saturn transfer')




