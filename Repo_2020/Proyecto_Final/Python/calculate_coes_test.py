
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts 

#time parameters
tspan=3600*10
dt=100.0

#central body
cb=pd.earth

if __name__ == '__main__':

	perts=null_perts()
	perts['J2']=True

	op=OP(t.tle2coes('HJ-2A.txt', degres=True),tspan,dt,coes=True,perts=perts, degres=True)

	op.calculate_coes()
	op.plot_coes(show_plot=True,hours=True,title='Orbita del satélite HJ-2A con perturbación J2')
	op.plot_3d(show_plot=True,title='Orbita del satélite HJ-2A con perturbación J2')