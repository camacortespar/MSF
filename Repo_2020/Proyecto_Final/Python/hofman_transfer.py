import math
import numpy as np
import datetime
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import tools as t
import planetary_data as pd
import spice_tools as st

cb=pd.earth

if __name__ == '__main__':
	
	#initial and final altitudes
	r0=0.0
	r1=100.0

	coes0=[6797.2+r0,0.0,51.6446,180.0,105.6661,51.6446]
	coes1=[6797.2+r1,0.0,51.6446,180.0,105.6661,51.6446]

	sc0,sc1,sc2,delta_v,t_transfer=t.hohmann_transfer(coes0=coes0,coes1=coes1,dt=1,altitude=False,propagate=True)

	t.plot_n_orbits([sc0.rs,sc1.rs,sc2.rs]
		,labels=['ISS orbita actual','ISS orbita final','ISS acelerando']
		,cb=cb,show_plot=True,title='Transferencia de la ISS a una órbita 100 km más alta')

	print('Delta V0=\t %.3f (km/s)' %delta_v[0])
	print('Delta Vf=\t %.3f (km/s)' %delta_v[1])
	print('time of transfer=\t %.1f seconds (%.2f hours)' %(t_transfer, t_transfer/3600))