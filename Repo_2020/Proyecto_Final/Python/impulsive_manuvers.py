import numpy as np
import datetime
import math as m
import matplotlib.pyplot as plt
from scipy.integrate import ode
import spiceypy as spice
from time import time as time
from multiprocessing import Pool
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
import tools as _t
import spice_tools as st
import spice_data as sd
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

#time parameters
tspan=3600*24*3
dt=100.0

#central body
cb=pd.earth

date0='2020-04-03'

if __name__ == '__main__':
	
	#calculate initial state 
	r0_apogee,v0_apogee=_t.coes2rv([cb['radius']+26600.0,0.74,35.0,180.0,0.0,0.0],degres=True)
	r0_perigee,v0_perigee=_t.coes2rv([cb['radius']+26600.0,0.74,35.0,0.0,0.0,0.0],degres=True)
	coes0_rotated=[cb['radius']+26600,0.74,35.0,0.0,45.0,0.0]

	#calculate circular velocity at initial position
	v_circ_apogee=(pd.earth['mu']/_t.norm(r0_apogee))**0.5
	v_circ_perigee=(pd.earth['mu']/_t.norm(r0_perigee))**0.5

	#calculate norm velocity vectors
	v_apogee_normed=_t.normed(v0_apogee)
	v_perigee_normed=_t.normed(v0_perigee)

	#calculate scape velocity
	esc_v_apogee=_t.esc_v(_t.norm(r0_apogee))
	esc_v_perigee=_t.esc_v(_t.norm(r0_perigee))

	#calculate current velocities
	v0_apogee_norm=_t.norm(v0_apogee)
	v0_perigee_norm=_t.norm(v0_perigee)

	#calculate scape trayectory vectors
	v0_apogee_scape=v_apogee_normed*esc_v_apogee
	v0_perigee_scape=v_perigee_normed*esc_v_perigee

	#print velocities and differences
	print('Current velocity at apogee:\t%.2f km/s' % v0_apogee_norm)
	print('Circular velocity at apogee:\t%.2f km/s'%v_circ_apogee)
	print('Scape velocity at v0_apogee:\t%.2f km/s'%esc_v_apogee)
	print('Delta V:\t%.2f km/s'%(esc_v_apogee-v0_apogee_norm))
	print()
	print('Current velocity at perigee:\t%.2f km/s' % v0_perigee_norm)
	print('Circular velocity at perigee:\t%.2f km/s'%v_circ_perigee)
	print('Scape velocity at perigee:\t%.2f km/s'%esc_v_perigee)
	print('Delta V:\t%.2f km/s'%(esc_v_perigee-v0_perigee_norm))
	print()
	
	#create new pert dictionary
	perts=null_perts()
	perts['n_bodies']=[pd.moon]

	#initial conditions
	state0_original=r0_apogee.tolist()+v0_apogee.tolist()
	state0_apogee=r0_apogee.tolist()+v0_apogee_scape.tolist()
	state0_perigee=r0_perigee.tolist()+v0_perigee_scape.tolist()

	#create instances of orbit propagator
	op_original=OP(state0_original,tspan,dt,perts=perts,date0=date0)
	op_apogee=OP(state0_apogee,tspan,dt,perts=perts,date0=date0)
	op_perigee=OP(state0_perigee,tspan,dt,perts=perts,date0=date0)
	op_rotated=OP(coes0_rotated,tspan,dt,coes=True,degres=True,perts=perts,date0=date0)

	#extract moon state vector
	moon_traj=op_apogee.perts['n_bodies'][0]['states'][:,:3]

	#plot all
	_t.plot_n_orbits([op_original.rs,op_apogee.rs,op_perigee.rs,op_rotated.rs,moon_traj],
		labels=['original','apoge','perigee','rotated','Moon'],show_plot=True)




