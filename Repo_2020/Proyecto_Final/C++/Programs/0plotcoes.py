import numpy as np
import plotingfunctions as t
if __name__ == '__main__':
	data0=np.loadtxt('OPcoes.dat',delimiter='	')
t.plot_coes(data0,hours=True,days=False,show_plot=True,save_plot=False,title='ISS COEs')
