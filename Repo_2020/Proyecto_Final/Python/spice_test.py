
from sys import path
path.append('')
import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
plt.style.use('dark_background')

import planetary_data as pd
import tools as _t
import spice_tools as st
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts


#central body
cb=pd.earth

#total steps of data
STEPS=100000

#reference frame of data
FRAME='ECLIPJ2000'

#observer of planetary bodies
OBSERVER='SUN'

if __name__ == '__main__':
	#load metakernel for solar system ephemerides
	spice.furnsh('spice_data/solar_system_kernel.mk')

	#get data from ephemerides file
	ids,names,tsc_sec,tsc_cal=st.get_objects('spice_data/de432s.bsp',display=True)

	#only include barycenters
	names=[f for f in names if 'BARYCENTER' in f]

	#create time array for ephemerides data for all bodies (assume all are the same)
	times=st.tc2array(tsc_sec[0],STEPS)

	#create empty list for all ephemerides data
	rs=[]

	#for each body in solar system
	for name in names:

		#add ephemerides data to list
		rs.append(st.get_ephemerides_data(name,times,FRAME,OBSERVER))

	_t.plot_n_orbits(rs,names,show_plot=True,cb=pd.sun)







