import numpy as np
import plotingfunctions as t
import planetary_data as pd
cb=pd.EARTH
if __name__ == '__main__':
	data0=np.loadtxt('OP.dat',delimiter='	')
	data0=data0[:,:3]
t.plot_n_orbits([data0,],labels=['ISS',],cb=pd.EARTH,show_plot=True,save_plot=False,title='Órbita del satélite HJ-2A')
