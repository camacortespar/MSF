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
import spice_data as sd

if __name__ == '__main__':

	spice.furnsh(sd.leap_seconds_kernel)

	sec2day=1/(3600*24)

	config={
	'planet0':'earth',
	'planet1':'venus',
	'departure0':'2005-01-01',
	'departure1':'2006-01-01',
	'arrival0':'2006-01-01',
	'arrival1':'2006-08-01',
	'mu':pd.sun['mu'],
	'step':1/sec2day,
	'frame':'ECLIPJ2000',
	'observer':'SOLAR SYSTEM BARYCENTER',
	'cuttof_v':20.0,
	'c3_levels':None,
	'vinf_levels':None,
	'tof_levels':None,
	'dv_levels':None,
	'dv_cmap':'RdPu_r',
	'figsize':(10,15),
	'lw':1.5,
	'title':'Porkchop plot of Jupiter to Saturn transfer',
	'title_dv':'Porkchop plot of Jupiter to Saturn transfer',
	'fontsize':15,
	'show':False,
	'filename':'venus2.png',
	'filename_dv':'venus3.png',
	'dpi':200,
	'load':False,
	}
	
	lt.interplanetary_porkchop_plot(config)
