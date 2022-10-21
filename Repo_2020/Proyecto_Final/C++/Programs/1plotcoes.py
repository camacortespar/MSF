import numpy as np
import plotingfunctions as t
if __name__ == '__main__':
	data1=np.loadtxt('OP1coes.dat',delimiter='	')
t.plot_coes(data1,hours=True,days=False,show_plot=True,save_plot=False,title='HJ-2A COEs')
