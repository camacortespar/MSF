import numpy as np
import plotingfunctions as t
if __name__ == '__main__':
	data2=np.loadtxt('OP2coes.dat',delimiter='	')
t.plot_coes(data2,hours=True,days=False,show_plot=True,save_plot=False,title='COSMOS2251 COEs')
